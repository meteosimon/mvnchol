load("../data/simdata.rda")
x  <- bamlssMVN:::log_dmvnchol_C(y, par)
xx <- bamlssMVN:::log_dmvnchol_ref(y, par)
all.equal(x, xx)

microbenchmark::microbenchmark(
  bamlssMVN:::log_dmvnchol_C(y, par),
  bamlssMVN:::log_dmvnchol_ref(y, par)
)

