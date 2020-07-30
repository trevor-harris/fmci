#' @export
fpc1 = function(h) {
  hfpc = refund::fpca.face(t(h), npc = 2)
  hs = t(hfpc$scores)
  hs[1,]
}
