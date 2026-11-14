// Rotation-shaped GEMM (U' G: U N x N-1, G N x S) with Eigen vs an external BLAS dgemm_ (linked at build time).
// Usage: ./gemm_<blas> N S reps
#include <Eigen/Dense>
#include <chrono>
#include <cstdio>
#include <cstdlib>
extern "C" void dgemm_(const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*,
                       const double*, const int*, const double*, double*, const int*);
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
int main(int argc, char** argv) {
  int N = atoi(argv[1]), S = atoi(argv[2]), reps = atoi(argv[3]); int M = N - 1;
  Eigen::MatrixXd U = Eigen::MatrixXd::Random(N, M), G = Eigen::MatrixXd::Random(N, S), R1(M, S), R2(M, S);
  double te = 1e9, tb = 1e9;
  for (int r = 0; r < reps; ++r) {
    double t0 = now(); R1.noalias() = U.transpose() * G; te = std::min(te, now() - t0);
    const double one = 1, zero = 0; t0 = now();
    dgemm_("T", "N", &M, &S, &N, &one, U.data(), &N, G.data(), &N, &zero, R2.data(), &M); tb = std::min(tb, now() - t0);
  }
  double gf = 2.0 * N * M * S / 1e9;
  printf("N %d S %d  eigen %.4f s (%.1f GF/s)  blas %.4f s (%.1f GF/s)  ratio eigen/blas %.2f  maxdiff %.2e\n", N, S, te, gf / te, tb, gf / tb, te / tb,
         (R1 - R2).cwiseAbs().maxCoeff() / R1.cwiseAbs().maxCoeff());
}
