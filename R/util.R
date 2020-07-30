#' @export
l1fit = function(y, cp) {
  cp = c(1, cp, length(y))

  yhat = y
  for(i in 2:length(cp)) {
    yhat[cp[i-1]:cp[i]] = mean(y[cp[i-1]:cp[i]])
  }
  yhat
}

#' @export
l1bic = function(y, cp) {
  yhat = l1fit(y, cp);
  n = length(y)
  n*log(var(y - yhat)) + sum(diff(yhat) != 0)*log(n);
}
