#' @export
proj_PELT = function(h) {
  round(poolC(sort(unique(c(changepoint::cpt.meanvar(ftv(h), method = "PELT")@cpts,
                            changepoint::cpt.meanvar(fpc1(h), method = "PELT")@cpts)))), 0)
}

#' @export
proj_FPOP = function(h) {
  round(poolC(sort(unique(c(robseg::Rob_seg.std(ftv(h), "L1", lambda = 5*log(ncol(h)))$t.est,
                            robseg::Rob_seg.std(fpc1(h), "L1", lambda = 5*log(ncol(h)))$t.est)))), 0)

}

#' @export
proj_WBS = function(h) {
  c(round(poolC(sort(unique(c(wbs::wbs(ftv(h))$cpt$cpt.ic$ssic.penalty,
                              wbs::wbs(fmean(h))$cpt$cpt.ic$ssic.penalty,
                              ncol(h))))), 0))
}
