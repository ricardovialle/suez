# Rscript cmp.R ref.rds new.rds: max |diff| per numeric column, row order/keys, NA and df identity
a = commandArgs(TRUE); A = as.data.frame(readRDS(a[1])); B = as.data.frame(readRDS(a[2]))
same_keys = nrow(A) == nrow(B) && identical(A$gene, B$gene) && identical(A$cis_snp, B$cis_snp)
num = names(A)[sapply(A, is.numeric)]
d = sapply(num, function(k) { x = abs(A[[k]] - B[[k]]); if (all(is.na(x))) 0 else max(x, na.rm = TRUE) })
na = sum(sapply(num, function(k) sum(is.na(A[[k]]) != is.na(B[[k]]))))
lcol = grep("^l", num, value = TRUE); bcol = grep("boot", lcol, value = TRUE); dcol = setdiff(lcol, bcol)
ok = same_keys && na == 0 && max(d[dcol]) <= 1e-9 && (!length(bcol) || max(d[bcol]) <= 1e-5) &&
  (!"df_interact" %in% num || identical(A$df_interact, B$df_interact))
cat(sprintf("rows %d keys %s NA-mismatch %d | max dl(data) %.2g dl(boot) %.2g | dbeta %.2g dse %.2g | identical %s => %s\n",
  nrow(B), same_keys, na, max(d[dcol]), if (length(bcol)) max(d[bcol]) else 0,
  if ("beta_geno" %in% num) d[["beta_geno"]] else NA, if ("se_geno" %in% num) d[["se_geno"]] else NA,
  identical(A, B), if (ok) "PASS" else "FAIL"))
