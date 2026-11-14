# Analytic eigen-decomposition of the default kernel K = Z Z' (Z = samples x individuals indicator) on the
# N-1 contrasts orthogonal to the intercept, as a SPARSE orthonormal basis (CSC lists, no Matrix package).
# Replaces kernel_contrast_eigen(tcrossprod(Z)) (O(N^3), N^2 memory) by O(N log N) work and memory.
#   null (0):      within-individual contrasts, hierarchical (Haar-like) splits of each individual's samples
#   size s:        between individuals of the same sample count s, Haar splits over those individuals
#   between sizes: k-1 vectors from the k x k problem on the size-class indicators (k = distinct sizes)
# Eigenvalues are exact (0, s) except the k-1 from the small k x k eigen-problem.

# orthonormal basis of 1^perp in R^m by recursive halving: list of (idx, val) vectors, m-1 of them, nnz ~ m log2 m
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

# dense version (tests / dense engine)
sp_dense = function(U) {
  M = matrix(0, U$dim[1], U$dim[2]); j = rep(seq_len(U$dim[2]), diff(U$p)); M[cbind(U$i + 1, j)] = U$x; M
}
