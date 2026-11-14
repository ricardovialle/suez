// C = A' B with A (k x m), B (k x n): Eigen vs external BLAS dgemm_. Usage: ./shape_<blas> m n k reps
#include <Eigen/Dense>
#include <chrono>
#include <cstdio>
#include <cstdlib>
extern "C" void dgemm_(const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*,
                       const double*, const int*, const double*, double*, const int*);
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
int main(int argc, char** argv) {
  int m = atoi(argv[1]), n = atoi(argv[2]), k = atoi(argv[3]), reps = atoi(argv[4]);
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(k, m), B = Eigen::MatrixXd::Random(k, n), C1(m, n), C2(m, n);
  double te = 1e9, tb = 1e9;
  for (int r = 0; r < reps; ++r) {
    double t0 = now(); C1.noalias() = A.transpose() * B; te = std::min(te, now() - t0);
    const double one = 1, zero = 0; t0 = now();
    dgemm_("T", "N", &m, &n, &k, &one, A.data(), &k, B.data(), &k, &zero, C2.data(), &m); tb = std::min(tb, now() - t0);
  }
  double gf = 2.0 * m * n * k / 1e9;
  printf("m %d n %d k %d  eigen %.2f GF/s  blas %.2f GF/s  eigen/blas time %.2f\n", m, n, k, gf / te, gf / tb, te / tb);
}
