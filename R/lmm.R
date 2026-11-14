## Fitting helpers for map_interaction_qtl: the model y = X b + e,
## e ~ N(0, sigma2_k * K + sigma2 * I), fitted by maximum likelihood
## (src/lmm.cpp, all SNPs of a gene per call) on data rotated by K's eigenvectors.
##
## The likelihood is that of the N-1 contrasts orthogonal to the intercept
## (REML for the intercept only). Every model has an intercept, so likelihood-
## ratio tests of the other effects are unchanged. This matters for a kernel
## with K 1 = 0, e.g. a GRM built from centred genotypes (plink, GCTA): on the
## full data the intercept fits that direction exactly and the likelihood grows
## without bound as sigma2 -> 0, so the reported log-likelihoods are arbitrary.

# Eigen-decomposition of K restricted to the contrasts: vectors = Q V
# (N x N-1, Q an orthonormal basis of the complement of 1). Refuses a kernel
# that is not positive semi-definite. Q is the Householder reflection
# H = I - 2 u u' of qr(1) without its first column, so Q'KQ = (HKH)[-1, -1] and
# QV = H [0; V] cost O(N^2) instead of dense products. The symmetric eigen-
# decomposition (eigen_solver()) is Eigen's (src/lmm.cpp), 2x faster than R's
# reference LAPACK at N = 1000-2000, or LAPACK through eigen() when R links an
# optimised BLAS/LAPACK: with OpenBLAS (1 thread) 4x faster than Eigen's at
# N = 2000 and 6.6x at N = 5000 (bench/proto/round11). Inside tied eigenvalues
# the two give different bases: data fits are identical, bootstrap draws for a
# given seed differ (same distribution).
kernel_contrast_eigen = function(K) {
  K = as.matrix(K)
  stopifnot(nrow(K) == ncol(K), is_symmetric(K, tol = 1e-8))
  n = nrow(K)
  u = -qr.Q(qr(matrix(1, n, 1)))[, 1]; u[1] = u[1] + 1; u = u / sqrt(sum(u^2))
  K = K - 2 * outer(drop(K %*% u), u)          # K H (R reuses outer()'s storage: one new N x N)
  # H K H, the eigen-decomposition of its contrast block and H [0; V] in one call, with R's arithmetic (u'(K H)
  # here, crossprod(u, V) by the same BLAS call inside): holds K + K H + U instead of ~7 N x N copies
  e = .Call(suez_contrast_eigen, K, u, drop(crossprod(u, K)), eigen_solver() == "lapack")
  rm(K)
  if (min(e$values) < -1e-6 * max(abs(e$values)))
    stop("sample_kernel is not positive semi-definite (min eigenvalue ",
         signif(min(e$values), 3), ").", call. = FALSE)
  e$values = pmax(e$values, 0)
  e
}

# options(suez.eigen_solver = "auto" | "eigen" | "lapack"). auto: LAPACK when R's BLAS or LAPACK (symlinks
# resolved, e.g. Debian alternatives) is an optimised library, else Eigen's solver (R's reference LAPACK).
eigen_solver = function() {
  s = match.arg(getOption("suez.eigen_solver", "auto"), c("auto", "eigen", "lapack"))
  if (s != "auto") return(s)
  libs = c(La_library(), tryCatch(extSoftVersion()[["BLAS"]], error = function(e) ""))
  libs = suppressWarnings(normalizePath(libs[nzchar(libs)], mustWork = FALSE))
  if (any(grepl("openblas|mkl|accelerate|veclib|blis|flexiblas|atlas", libs, ignore.case = TRUE))) "lapack" else "eigen"
}

# isSymmetric(unname(K), tol) (all.equal: mean relative difference over the unequal entries), one column
# at a time: isSymmetric holds ~4 N x N temporaries.
is_symmetric = function(K, tol) {
  d = a = m = 0
  for (j in seq_len(ncol(K))) {
    x = K[, j]; y = K[j, ]; u = x != y
    if (anyNA(u)) return(FALSE)
    d = d + sum(abs(x[u] - y[u])); a = a + sum(abs(x[u])); m = m + sum(u)
  }
  m == 0 || (if (a / m > tol) d / a else d / m) <= tol
}

# Null directions of the kernel (eigenvalue ~ 0, e.g. contrasts between samples
# of the same individual with the default kernel) only have variance sigma2. If
# the fixed effects fit all of them exactly the likelihood is unbounded as
# sigma2 -> 0. `need` = columns the alternative models add that can vary within
# those directions.
check_identifiable = function(eig, X0, need = 1) {
  null = eig$values <= 1e-8 * max(eig$values)
  m = sum(null)
  if (m == 0) return(invisible(Inf))
  d = if (ncol(X0)) svd(crossprod_U(eig$vectors, X0)[null, , drop = FALSE], 0, 0)$d else 0
  r = sum(d > 1e-8 * max(d, 1))   # svd rank: qr's rank would count a ~0 column
  if (m - r < need) stop(
    "the sample kernel leaves ", m, " direction(s) with no kernel variance (e.g. contrasts ",
    "between samples of the same individual); the condition design fits ", r, " and the ",
    "tested effects need ", need, " more, so sigma2 is not identifiable.", call. = FALSE)
  invisible(m - r)
}

# Centre and scale columns (the likelihood is invariant; L-BFGS stops short on
# raw-scale columns such as age).
standardize_cols = function(M) {
  if (ncol(M)) {
    M = scale(M)
    attr(M, "scaled:center") = attr(M, "scaled:scale") = NULL
  }
  M
}
haar = function(m) {
  out = list()
  rec = function(lo, hi) {   # items lo..hi
    m = hi - lo + 1; if (m < 2) return()
    a = m %/% 2; L = lo:(lo + a - 1); R = (lo + a):hi; nl = length(L); nr = length(R)
    out[[length(out) + 1]] <<- list(idx = c(L, R), val = c(rep(sqrt(nr / (m * nl)), nl), rep(-sqrt(nl / (m * nr)), nr)))
    rec(lo, lo + a - 1); rec(lo + a, hi)
  }
  rec(1, m); out
}

zz_contrast_basis = function(individual) {
  ind = match(individual, unique(individual)); N = length(ind); J = max(ind)
  nj = tabulate(ind, J); samples = split(seq_len(N), ind)          # samples of each individual, increasing
  ti = list(); tx = list(); lam = numeric(0); col = 0
  add = function(i, x, l) { col <<- col + 1; ti[[col]] <<- i; tx[[col]] <<- x; lam[col] <<- l }
  hcache = list()
  H = function(m) { key = as.character(m); if (is.null(hcache[[key]])) hcache[[key]] <<- haar(m); hcache[[key]] }
  # 1. within-individual (eigenvalue 0)
  for (j in which(nj >= 2)) { sj = samples[[j]]; for (v in H(nj[j])) add(sj[v$idx], v$val, 0) }
  # 2. between individuals of one size class (eigenvalue s): x = Z a / sqrt(s)
  sizes = sort(unique(nj))
  for (s in sizes) {
    js = which(nj == s); if (length(js) < 2) next
    for (v in H(length(js))) {
      ii = unlist(samples[js[v$idx]], use.names = FALSE)
      add(ii, rep(v$val / sqrt(s), each = s), s)
    }
  }
  # 3. between size classes: e_s = class indicator / sqrt(N_s); P diag(s) P on q^perp, q = sqrt(N_s)
  k = length(sizes)
  if (k >= 2) {
    Ns = sapply(sizes, function(s) s * sum(nj == s)); q = sqrt(Ns)
    Qk = qr.Q(qr(cbind(q / sqrt(sum(Ns)), diag(k))))[, -1, drop = FALSE]   # orthonormal basis of q^perp
    e = eigen(crossprod(Qk, sizes * Qk), symmetric = TRUE)
    W = Qk %*% e$vectors
    cls = match(nj[ind], sizes)
    for (r in seq_len(k - 1)) add(seq_len(N), W[cls, r] / sqrt(Ns[cls]), e$values[r])
  }
  stopifnot(col == N - 1)
  o = order(lam)                                  # ascending, as the base solver; ties in construction order
  ti = ti[o]; tx = tx[o]; lam = lam[o]
  len = lengths(ti); p = c(0L, cumsum(len))
  # within a column the rows must be increasing for CSC
  i = integer(sum(len)); x = numeric(sum(len))
  for (c in seq_along(ti)) { r = order(ti[[c]]); sl = (p[c] + 1):p[c + 1]; i[sl] = ti[[c]][r] - 1L; x[sl] = tx[[c]][r] }
  list(U = list(i = as.integer(i), p = as.integer(p), x = x, dim = c(N, N - 1L)), values = lam)
}

# U'M for a dense U or the sparse structured basis
crossprod_U = function(U, M) if (is.list(U)) .Call(suez_sp_crossprod, U, as.matrix(M)) else crossprod(U, M)

# C++ threads per gene fit (per doMC worker): options(suez.threads = n), or the environment variable SUEZ_THREADS;
# default 1. Results are identical() at any thread count (SNP blocks are the parallel tasks; bootstrap normals are
# drawn on the R thread in the same order). Useful for few genes at large N, or to use cores with less RAM than
# forked workers (e.g. registerDoMC(2) with options(suez.threads = 4)).
suez_threads = function() {
  t = suppressWarnings(as.integer(getOption("suez.threads", Sys.getenv("SUEZ_THREADS", "1"))))
  if (length(t) != 1 || is.na(t) || t < 1) stop("options(suez.threads) must be a positive integer.", call. = FALSE)
  t
}
