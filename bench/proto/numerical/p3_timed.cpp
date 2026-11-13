// [[Rcpp::depends(RcppEigen)]]
// Exact ML fits for map_interaction_qtl, all cis-SNPs of one gene per call.
// Model on rotated data: Uy = UX b + e, e ~ N(0, diag(sigma2_k * lambda + sigma2)). Writing
// sigma2_k = s h, sigma2 = s (1 - h), b and s have closed forms for each h in [0, 1], so each fit is a 1-D
// maximisation over h. The search follows GEMMA / FaST-LMM: a grid uniform in x = log10(delta),
// delta = sigma2_k / sigma2 = h / (1 - h), evaluated for every SNP of the gene together (shared condition
// terms once, genotype terms as matrix products across SNPs); then safeguarded Newton steps with analytic
// derivatives around every local maximum of the grid (the likelihood can be bimodal), keeping the best.
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <chrono>
static double TT[8];
static double now(){return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();}
// [[Rcpp::export]]
Rcpp::NumericVector timers(bool reset){Rcpp::NumericVector v(TT,TT+8); if(reset) for(auto&x:TT)x=0; return v;}
using Eigen::ArrayXd; using Eigen::ArrayXXd; using Eigen::MatrixXd; using Eigen::VectorXd; using Eigen::VectorXi;
static const double NEG_INF = -std::numeric_limits<double>::infinity();
static const double L10 = 2.302585092994046;   // log(10)

// h grid: h = 0, then x = log10(delta) uniform on [-5, 5] (80 points, spacing 0.13), then h = 1. Uniform in
// x, the grid suits kernels of any scale; in h a fixed grid misses narrow peaks when the kernel's
// eigenvalues span orders of magnitude (e.g. a learned cov_suez). Peaks beyond |x| = 5 are reached from
// the end points (Newton works up to |x| = 12).
static std::vector<double> make_grid() {
  const int nx = 80; const double xl = 5;
  std::vector<double> h(1, 0.0);
  for (int i = 0; i < nx; ++i) h.push_back(1 / (1 + std::pow(10.0, -(-xl + 2 * xl * i / (nx - 1.0)))));
  h.push_back(1); return h;
}
static const std::vector<double> HS = make_grid();
static const int NG = HS.size();
static const int MAX_PEAKS = 3;
static long FALLBACKS = 0;
// [[Rcpp::export]]
double fallbacks() { double f = FALLBACKS; FALLBACKS = 0; return f; }   // local maxima of the grid refined per fit


// p3: eigenvalue groups. Contrasts with tied eigenvalues (|lam_t - lam_u| within 1e-12 max|lam|: the default
// kernel has 2 distinct values, a learned cov_suez ~20) share the weight w_t(h), so every weighted sum the
// likelihood needs, sum_t w_t u_t v_t, equals sum_g w_g S_g(u, v) with the group moments S_g(u, v) =
// sum_{t in g} u_t v_t, computed once per gene / SNP. The grid (all h at once: Wg (NG x G) * moments) and every
// Newton step then cost O(G) instead of O(n). Rows are permuted so groups are contiguous; the bootstrap
// draws stay in the original contrast order.
struct Groups {
  int n, G; ArrayXd lam, cnt, c; std::vector<int> perm, start; double lmin, lmax; MatrixXd W; VectorXd logdet; std::vector<char> bad;
  explicit Groups(const ArrayXd& l) : n(l.size()) {
    perm.resize(n); for (int i = 0; i < n; ++i) perm[i] = i;
    std::stable_sort(perm.begin(), perm.end(), [&](int a, int b) { return l[a] < l[b]; });
    lmin = l.minCoeff(); lmax = l.maxCoeff(); const double tol = 1e-12 * std::max(std::fabs(lmin), std::fabs(lmax));
    std::vector<double> lg, cn;
    for (int i = 0; i < n; ++i) {
      if (i == 0 || l[perm[i]] - l[perm[start.back()]] > tol) { start.push_back(i); lg.push_back(0); cn.push_back(0); }
      lg.back() += l[perm[i]]; cn.back() += 1;
    }
    G = start.size(); start.push_back(n);
    lam.resize(G); cnt.resize(G); for (int g = 0; g < G; ++g) { cnt[g] = cn[g]; lam[g] = lg[g] / cn[g]; }
    c = lam - 1; W.resize(NG, G); logdet.resize(NG); bad.resize(NG);
    for (int i = 0; i < NG; ++i) {
      bad[i] = !(HS[i] * lmin + (1 - HS[i]) > 0 && HS[i] * lmax + (1 - HS[i]) > 0);   // = any(d <= 0): d is linear in lam
      ArrayXd d = HS[i] * lam + (1 - HS[i]);
      if (bad[i]) { W.row(i).setZero(); logdet[i] = NAN; continue; }
      W.row(i) = d.inverse().matrix().transpose(); logdet[i] = (cnt * d.log()).sum();
    }
  }
  bool is_bad(double h) const { return !(h * lmin + (1 - h) > 0 && h * lmax + (1 - h) > 0); }
  MatrixXd permute(const MatrixXd& M) const { MatrixXd R(M.rows(), M.cols()); for (int i = 0; i < n; ++i) R.row(i) = M.row(perm[i]); return R; }
  // group sums of u .* v
  void moments(const double* u, const double* v, double* out) const {
    for (int g = 0; g < G; ++g) { double acc = 0; for (int t = start[g]; t < start[g + 1]; ++t) acc += u[t] * v[t]; out[g] = acc; }
  }
};

// pair index among N columns (a <= b), row-major upper triangle
static inline int pidx(int a, int b, int N) { if (a > b) std::swap(a, b); return a * N - a * (a - 1) / 2 + (b - a); }

// Group moments of the columns [X0, B_1..B_k, Y] for every SNP: Mom (G x npair*S), SNP-major; the X0 x X0
// pairs (and X0 x Y, Y x Y with a single y) are the same for every SNP and computed once.
struct Moments {
  int p0, k, N, iy, np, S; bool shared_y; MatrixXd M;
  Moments(const Groups& gr, const MatrixXd& X0, const std::vector<MatrixXd>& B, const MatrixXd& Y, int S_)
      : p0(X0.cols()), k(B.size()), N(p0 + k + 1), iy(p0 + k), np(N * (N + 1) / 2), S(S_), shared_y(Y.cols() == 1), M(gr.G, np * S_) {
    auto col = [&](int a, int s) { return a < p0 ? X0.col(a).data() : a < iy ? B[a - p0].col(s).data() : Y.col(shared_y ? 0 : s).data(); };
    auto shared = [&](int a) { return a < p0 || (a == iy && shared_y); };
    VectorXd tmp(gr.G);
    for (int a = 0; a < N; ++a) for (int b = a; b < N; ++b) {
      const int t = pidx(a, b, N);
      if (shared(a) && shared(b)) { gr.moments(col(a, 0), col(b, 0), tmp.data()); for (int s = 0; s < S; ++s) M.col(s * np + t) = tmp; }
      else for (int s = 0; s < S; ++s) gr.moments(col(a, s), col(b, s), M.col(s * np + t).data());
    }
  }
  // moments of the model [cols..., Y] of SNP s, as G x (p+1)(p+2)/2 in the same pair layout
  MatrixXd model(int s, const std::vector<int>& cols) const {
    const int p = cols.size(), Nm = p + 1; MatrixXd R(M.rows(), Nm * (Nm + 1) / 2);
    auto c = [&](int a) { return a < p ? cols[a] : iy; };
    for (int a = 0; a < Nm; ++a) for (int b = a; b < Nm; ++b) R.col(pidx(a, b, Nm)) = M.col(s * np + pidx(c(a), c(b), N));
    return R;
  }
};

struct Fit { double value = NAN, s2k = NAN, s2 = NAN; VectorXd b; bool ok = false; };

// Profile log-likelihood and derivatives as in p0, from group moments MM (G x pairs of [X, y]):
//   A = sum_g w_g Sxx_g, rhs = sum_g w_g Sxy_g, rr_g = |r|^2 on group g, g = sum w_g rr_g,
//   Q1 = sum c_g w_g^2 rr_g, Q2 = sum c_g^2 w_g^3 rr_g - Xv' A^-1 Xv, Xv = sum c_g w_g^2 (Sxy_g - Sxx_g b).
struct Profile {
  const Groups& gr; int n; ArrayXd d, w, rr; MatrixXd A; VectorXd b, rhs, mw, coef, Xv, tmp; Eigen::LDLT<MatrixXd> L; double s = NAN;
  explicit Profile(const Groups& g) : gr(g), n(g.n) {}
  double operator()(const MatrixXd& MM, double h, double* f1 = nullptr, double* f2 = nullptr) {
    const int Nm = (int)std::lround((std::sqrt(8.0 * MM.cols() + 1) - 1) / 2), p = Nm - 1;
    if (gr.is_bad(h)) { if (f1) *f1 = *f2 = NAN; return NEG_INF; }
    d = h * gr.lam + (1 - h); w = d.inverse();
    mw.noalias() = MM.transpose() * w.matrix();
    coef.setZero(MM.cols()); coef[pidx(p, p, Nm)] = 1;
    if (p) {
      A.resize(p, p); rhs.resize(p);
      for (int a = 0; a < p; ++a) { rhs[a] = mw[pidx(a, p, Nm)]; for (int e = 0; e <= a; ++e) A(a, e) = A(e, a) = mw[pidx(e, a, Nm)]; }
      L.compute(A); b = L.solve(rhs);
      for (int a = 0; a < p; ++a) { coef[pidx(a, p, Nm)] = -2 * b[a]; for (int e = a; e < p; ++e) coef[pidx(a, e, Nm)] = (a == e ? 1 : 2) * b[a] * b[e]; }
    } else b.resize(0);
    rr = (MM * coef).array();
    double g = (w * rr).sum(); s = g / n;
    if (!(s > 0)) { if (f1) *f1 = *f2 = NAN; return std::numeric_limits<double>::infinity(); }
    if (f1) {
      ArrayXd cw2 = gr.c * w.square();
      double Q1 = (cw2 * rr).sum(), Q2 = (cw2 * gr.c * w * rr).sum();
      if (p) {
        mw.noalias() = MM.transpose() * cw2.matrix(); Xv.resize(p);
        for (int a = 0; a < p; ++a) { double x = mw[pidx(a, p, Nm)]; for (int e = 0; e < p; ++e) x -= mw[pidx(std::min(a, e), std::max(a, e), Nm)] * b[e]; Xv[a] = x; }
        tmp = L.solve(Xv); Q2 -= Xv.dot(tmp);
      }
      *f1 = 0.5 * n * Q1 / g - 0.5 * (gr.cnt * gr.c * w).sum();
      *f2 = -0.5 * n * (2 * Q2 / g - Q1 * Q1 / (g * g)) + 0.5 * (gr.cnt * (gr.c * w).square()).sum();
    }
    return -0.5 * n - 0.5 * (n * std::log(s) + (gr.cnt * d.log()).sum());
  }
};
// Safeguarded Newton in x = log10(delta) on the bracket around grid point `best`: h = 1 / (1 + 10^-x),
// dh/dx = ln10 h (1 - h), d2h/dx2 = ln10 (1 - 2h) dh/dx. Steps leaving the bracket become bisections.
// Returns the best h found and its value.
static std::pair<double, double> newton(Profile& P, const MatrixXd& MM,
                                        const double* v, int best) {
  const double XCAP = 12;
  auto X_of = [&](double h) { return h <= 0 ? -XCAP : h >= 1 ? XCAP : std::max(-XCAP, std::min(XCAP, std::log10(h / (1 - h)))); };
  auto H = [](double x) { return 1 / (1 + std::pow(10.0, -x)); };
  double hb = HS[best], vb = v[best], f1, f2;
  double a = HS[std::max(0, best - 1)], c = HS[std::min(NG - 1, best + 1)];
  if (!std::isfinite(vb) || !(c > a)) return {hb, vb};
  double lo = X_of(a), hi = X_of(c), x = X_of(hb);
  for (int it = 0; it < 60; ++it) {
    double h = H(x), fv = P(MM, h, &f1, &f2);
    if (fv > vb) { vb = fv; hb = h; }
    double dh = L10 * h * (1 - h), d2h = L10 * (1 - 2 * h) * dh;
    double g1 = f1 * dh, g2 = f2 * dh * dh + f1 * d2h;
    if (!std::isfinite(g1)) break;
    if (g1 > 0) lo = std::max(lo, x); else hi = std::min(hi, x);       // the maximum is uphill
    double xn = g2 < 0 ? x - g1 / g2 : NAN;
    if (!std::isfinite(xn) || xn <= lo || xn >= hi) xn = 0.5 * (lo + hi);
    if (std::fabs(xn - x) < 1e-9 || hi - lo < 1e-10) {
      double hh = H(xn), fv2 = P(MM, hh);
      if (fv2 > vb) { vb = fv2; hb = hh; }
      break;
    }
    x = xn;
  }
  return {hb, vb};
}

// Refine every local maximum of the grid values v (up to MAX_PEAKS, highest first) and return the best fit.
static Fit refine(Profile& P, const MatrixXd& MM, const double* v, int best) {
  std::vector<int> peaks;
  for (int i = 0; i < NG; ++i)
    if (std::isfinite(v[i]) && (i == 0 || v[i] >= v[i - 1]) && (i == NG - 1 || v[i] >= v[i + 1])) peaks.push_back(i);
  std::sort(peaks.begin(), peaks.end(), [&](int p, int q) { return v[p] > v[q]; });
  if (peaks.empty()) peaks.push_back(best);
  double hb = HS[peaks[0]], vb = NEG_INF;
  for (size_t j = 0; j < peaks.size() && j < (size_t)MAX_PEAKS; ++j) {
    auto m = newton(P, MM, v, peaks[j]);
    if (m.second > vb) { vb = m.second; hb = m.first; }
  }
  Fit f; f.value = P(MM, hb); f.s2k = P.s * hb; f.s2 = P.s * (1 - hb); f.b = P.b;
  f.ok = std::isfinite(f.value);
  return f;
}

// Grid values (S x NG) for the nested models with the first kk blocks, kk in kks, as p2 (Schur complement on
// the X0 block, LDLT fallback on a near-singular pivot) with moments R = Wg * Mom.
static std::vector<MatrixXd> grid_values(const Groups& gr, const Moments& mo, const std::vector<int>& kks) {
  const int p0 = mo.p0, k = mo.k, N = mo.N, iy = mo.iy, np = mo.np, S = mo.S, n = gr.n;
  MatrixXd R = gr.W * mo.M;                                 // NG x np*S
  std::vector<MatrixXd> V(kks.size(), MatrixXd(S, NG));
  std::vector<MatrixXd> A0i(NG);
  for (int i = 0; i < NG; ++i) if (!gr.bad[i]) {
    MatrixXd A00(p0, p0); for (int a = 0; a < p0; ++a) for (int b = 0; b < p0; ++b) A00(a, b) = R(i, pidx(a, b, N));
    A0i[i] = A00.ldlt().solve(MatrixXd::Identity(p0, p0));
  }
  VectorXd r0(p0), z(p0), u(k), c(k); MatrixXd a(p0, k), T(p0, k), Sc(k, k);
  for (int i = 0; i < NG; ++i) {
    if (gr.bad[i]) { for (auto& v : V) v.col(i).setConstant(NEG_INF); continue; }
    for (int s = 0; s < S; ++s) {
      auto M = [&](int x, int y) { return R(i, s * np + pidx(x, y, N)); };
      double yy = M(iy, iy), q = 0;
      if (p0) { for (int e = 0; e < p0; ++e) r0[e] = M(e, iy); z.noalias() = A0i[i] * r0; q = r0.dot(z); }
      for (size_t mi = 0; mi < kks.size(); ++mi) {
        const int kk = kks[mi]; double rss = yy - q; bool fb = false;
        if (kk) {
          for (int j = 0; j < kk; ++j) {
            for (int e = 0; e < p0; ++e) a(e, j) = M(e, p0 + j);
            c[j] = M(p0 + j, iy);
            for (int l = 0; l <= j; ++l) Sc(j, l) = M(p0 + j, p0 + l);
          }
          auto aa = a.leftCols(kk); auto TT = T.leftCols(kk);
          if (p0) TT.noalias() = A0i[i] * aa;
          for (int j = 0; j < kk; ++j) {
            if (Sc(j, j) == 0) { for (int l = 0; l < j; ++l) Sc(j, l) = 0; Sc(j, j) = 1; u[j] = 0; continue; }   // absent column
            u[j] = c[j] - (p0 ? aa.col(j).dot(z) : 0.0);
            if (p0) for (int l = 0; l <= j; ++l) Sc(j, l) -= aa.col(j).dot(TT.col(l));
          }
          for (int j = 0; j < kk; ++j) {
            double dj = Sc(j, j); const double dj0 = dj;
            for (int l = 0; l < j; ++l) dj -= Sc(j, l) * Sc(j, l);
            if (!(dj > 1e-10 * std::max(dj0, 1e-300))) { fb = true; break; }
            Sc(j, j) = std::sqrt(dj);
            for (int r = j + 1; r < kk; ++r) { double x = Sc(r, j); for (int l = 0; l < j; ++l) x -= Sc(r, l) * Sc(j, l); Sc(r, j) = x / Sc(j, j); }
            double x = u[j]; for (int l = 0; l < j; ++l) x -= Sc(j, l) * u[l]; u[j] = x / Sc(j, j); rss -= u[j] * u[j];
          }
          if (fb) {
            const int p = p0 + kk; MatrixXd A(p, p); VectorXd rhs(p);
            for (int e = 0; e < p; ++e) { rhs[e] = M(e, iy); for (int f = 0; f < p; ++f) A(e, f) = M(e, f); }
            for (int j = 0; j < kk; ++j) if (A(p0 + j, p0 + j) == 0) { A(p0 + j, p0 + j) = 1; rhs[p0 + j] = 0; }
            rss = yy - rhs.dot(A.ldlt().solve(rhs)); ++FALLBACKS;
          }
        }
        double sv = rss / n;
        V[mi](s, i) = sv > 0 ? -0.5 * n - 0.5 * (n * std::log(sv) + gr.logdet[i]) : std::numeric_limits<double>::infinity();
      }
    }
  }
  return V;
}

static int argmax(const double* v) { int b = 0; for (int i = 1; i < NG; ++i) if (v[i] > v[b]) b = i; return b; }

static bool standardize(Eigen::Ref<VectorXd> x, double& sd) {
  const int n = x.size(); x.array() -= x.mean(); sd = std::sqrt(x.squaredNorm() / (n - 1));
  if (!(sd > 0)) return false;
  x /= sd; return true;
}

static Rcpp::List fit_gene(const Eigen::Map<MatrixXd>& U, const VectorXd& lam_, const VectorXd& Uy_, const MatrixXd& UX0_,
                           const Eigen::Map<MatrixXd>& G, const MatrixXd& C, bool boot) {
  const int S = G.cols(), q = C.cols(), p0 = UX0_.cols(), n = Uy_.size();
  const ArrayXd lam = lam_.array();
  Groups gr(lam); Profile PR(gr);
  const MatrixXd UX0 = gr.permute(UX0_); MatrixXd Y1 = gr.permute(MatrixXd(Uy_));
  std::vector<int> c0(p0); for (int e = 0; e < p0; ++e) c0[e] = e;
  auto run = [&](const Moments& mo, int s, const std::vector<int>& cols, const double* v) {
    MatrixXd MM = mo.model(s, cols); return std::make_pair(refine(PR, MM, v, argmax(v)), MM);
  };

  Moments m0(gr, UX0, {}, Y1, 1);
  MatrixXd V0 = grid_values(gr, m0, {0})[0];
  Fit f0 = run(m0, 0, c0, V0.data()).first;
  if (!f0.ok) Rcpp::stop("no-genotype model did not converge");

  double t0=now();
  MatrixXd Gs = G; VectorXd gsd(S);
  for (int s = 0; s < S; ++s) standardize(Gs.col(s), gsd[s]);
  std::vector<MatrixXd> blocks(1, gr.permute(MatrixXd(U.transpose() * Gs)));
  std::vector<std::vector<bool>> iok(q, std::vector<bool>(S)); VectorXi df = VectorXi::Zero(S);
  for (int j = 0; j < q; ++j) {
    MatrixXd I = G.array().colwise() * C.col(j).array();
    for (int s = 0; s < S; ++s) { double sd; iok[j][s] = standardize(I.col(s), sd); if (!iok[j][s]) I.col(s).setZero(); else ++df[s]; }
    blocks.push_back(gr.permute(MatrixXd(U.transpose() * I)));
  }
  TT[0]+=now()-t0; t0=now();
  const MatrixXd& UG = blocks[0];
  const int K = blocks.size();
  auto cols1 = [&](int) { std::vector<int> c = c0; c.push_back(p0); return c; };
  auto cols2 = [&](int s) { std::vector<int> c = cols1(s); for (int j = 0; j < q; ++j) if (iok[j][s]) c.push_back(p0 + 1 + j); return c; };
  auto ok_nested = [](Fit& f, const Fit& g0) { if (!f.ok || f.value < g0.value - 1e-6) f.ok = false; };

  Rcpp::NumericMatrix out(S, 7);
  std::fill(out.begin(), out.end(), NA_REAL);
  Moments md(gr, UX0, blocks, Y1, S);
  std::vector<MatrixXd> VV = grid_values(gr, md, q ? std::vector<int>{1, K} : std::vector<int>{1});
  std::vector<Fit> f1(S); TT[1]+=now()-t0; t0=now();
  for (int s = 0; s < S; ++s) {
    out(s, 2) = df[s];
    VectorXd v1 = VV[0].row(s).transpose(); f1[s] = run(md, s, cols1(s), v1.data()).first; ok_nested(f1[s], f0);
    if (!f1[s].ok) continue;
    out(s, 0) = f1[s].value;
    out(s, 3) = f1[s].b[p0] / gsd[s]; out(s, 4) = std::sqrt(PR.s * PR.L.solve(VectorXd::Unit(p0 + 1, p0))[p0]) / gsd[s];   // Vb = s A^-1
    if (df[s] > 0) {
      VectorXd v2 = VV[1].row(s).transpose(); Fit f2 = run(md, s, cols2(s), v2.data()).first; ok_nested(f2, f1[s]);
      if (f2.ok) out(s, 1) = f2.value;
    }
  }
  TT[2]+=now()-t0; t0=now();
  if (!boot) return Rcpp::List::create(Rcpp::Named("l0") = f0.value, Rcpp::Named("res") = out);

  // bootstrap: drawn in the original contrast order (SNP by SNP), stored at the permuted rows
  Rcpp::RNGScope rng;
  std::vector<int> pos(n); for (int i = 0; i < n; ++i) pos[gr.perm[i]] = i;
  MatrixXd Yb = MatrixXd::Zero(n, S);
  for (int s = 0; s < S; ++s) {
    if (!f1[s].ok) continue;
    for (int i = 0; i < n; ++i) Yb(pos[i], s) = std::sqrt(f1[s].s2k * lam[i] + f1[s].s2) * R::norm_rand();
    Yb.col(s) += UX0 * f1[s].b.head(p0) + UG.col(s) * f1[s].b[p0];
  }
  TT[3]+=now()-t0; t0=now();
  Moments mb(gr, UX0, blocks, Yb, S); TT[4]+=now()-t0; t0=now();
  std::vector<MatrixXd> WW = grid_values(gr, mb, q ? std::vector<int>{0, 1, K} : std::vector<int>{0, 1});
  TT[5]+=now()-t0; t0=now();
  for (int s = 0; s < S; ++s) {
    if (!f1[s].ok) continue;
    VectorXd v = WW[0].row(s).transpose(); Fit g0 = run(mb, s, c0, v.data()).first;
    if (!g0.ok) continue;
    v = WW[1].row(s).transpose(); Fit g1 = run(mb, s, cols1(s), v.data()).first; ok_nested(g1, g0);
    if (!g1.ok) continue;
    out(s, 5) = g1.value;
    if (df[s] > 0) {
      v = WW[2].row(s).transpose(); Fit g2 = run(mb, s, cols2(s), v.data()).first; ok_nested(g2, g1);
      if (g2.ok) out(s, 6) = g2.value;
    }
  }
  TT[6]+=now()-t0;
  return Rcpp::List::create(Rcpp::Named("l0") = f0.value, Rcpp::Named("res") = out);
}

// [[Rcpp::export]]
Rcpp::List fit_gene_proto(Eigen::Map<Eigen::MatrixXd> U, Eigen::VectorXd lambda, Eigen::VectorXd Uy, Eigen::MatrixXd UX0,
                          Eigen::Map<Eigen::MatrixXd> G, Eigen::MatrixXd C, bool boot) { return fit_gene(U, lambda, Uy, UX0, G, C, boot); }
