#' @export
ks_cdf = function(x, n = 20) {
  if(x < 0.05) return(0)
  1 - 2*(sum(sapply(1:n, function(k) ((-1)^(k-1)) * exp(-2*(k^2)*(x^2)))))
}

#' @export
cplink = function(yhat, eps = sqrt(length(yhat)), cp = NULL) {
  # eps = eps*sqrt(length(yhat))

  if(any(is.null(cp))) {
    st = diff(yhat)
    cp = which(st != 0)
  }

  cpb = list()
  ncp = length(cp)

  i = 1
  while(i <= ncp) {
    group = cp[i]
    if(i != ncp) {

      # link all CP within reach to form the group
      for(j in (i+1):ncp) {
        if((cp[j] - cp[j-1]) < eps) {
          group = c(group, cp[j])
        } else {
          break
        }
      }
      cpb = append(cpb, list(group))

      # get the index number so the associated weight can be fetched
      ind = which(cp %in% group)

      i = i + length(ind)
    } else {
      i = i + 1
      cpb = append(cpb, list(group))
    }
  }
  append(1L, append(cpb, length(yhat)))
}

#' @export
cpfilt = function(y, yhat, cpb, alpha = 1e-3) {
  cps = c()
  pval = c()

  if(length(cpb) == 2) {
    # return(list(cps = cps, pval = pval))
    return(cps)
  }

  for(i in 2:(length(cpb)-1)) {
    rng = seq(cpb[[i-1]][length(cpb[[i-1]])], cpb[[i+1]][1])
    yrng = y[rng]

    # aux variables
    yres = yrng - yhat[rng]
    ylen = length(yrng)
    ymu = mean(yrng)
    ysd = sd(yres)

    # test statistic and pval
    tn = abs(1/sqrt(ylen) * cumsum(yrng - ymu) / ysd)
    pv = 1-ks_cdf(max(tn))

    cps = c(cps, rng[which.max(tn)])
    pval = c(pval, pv)
  }
  pval = p.adjust(pval, "BH")
  unique(sort(cps[pval <= alpha]))
  # pval = pval[pval <= alpha]
  # list(cps = cps[pval <= alpha], pval = pval[pval <= alpha])
}

#' @export
eps_opt = function(y, alpha = 1e-3, lam = 1, eps0 = sqrt(length(y))) {
  len = 100
  bics = rep(0, length.out = len)
  eps_seq = seq(0.1, 10, length.out = len)

  for(b in 1:length(bics)) {
    cs = mci(y, alpha = alpha, opt_level = 0, eps = eps_seq[b] * eps0, lam = lam)

    if(is.null(cs)) {
      break
    } else {
      bics[b] = l1bic(y, cs)
    }
  }
  eps_seq[which.min(bics)] * eps0
}

#' @export
lam_opt = function(y, alpha = 1e-3, eps = 1, lam0 = sqrt(length(y))) {
  len = 50
  bics = rep(0, length.out = len)
  c_seq = seq(0.2, 5, length.out = len)

  for(b in 1:len) {
    cs = mci(y, alpha = alpha, opt_level = 0, eps0 = eps, lam0 = c_seq[b] * lam0)

    if(is.null(cs)) {
      break
    } else {
      bics[b] = l1bic(y, cs)
    }
  }
  c_seq[which.min(bics)] * lam0
}

#' @export
param_opt = function(y, alpha = 1e-3, eps0 = sqrt(length(y)), lam0 = sqrt(length(y))) {
  len = 25
  bics = matrix(0, len, len)

  c_seq = seq(0.01, 2, length.out = len)
  k_seq = seq(0.01, 2, length.out = len)

  for(c in 1:len) {
    for(k in 1:len) {
      cs = mci(y, alpha = alpha, opt_level = 0, eps0 = c_seq[c] * eps0, lam0 = k_seq[k] * lam0)
      bics[c, k] = l1bic(y, cs)
    }
  }

  ck = which(bics == min(bics), arr.ind = T)[1,]
  eps_hat = c_seq[ck[1]] * eps0
  lam_hat = k_seq[ck[2]] * lam0
  return(c(eps_hat, lam_hat))
}

#' @export
mci = function(y, alpha = 1e-3, opt_level = 1, lam0 = sqrt(length(y)), eps0 = sqrt(length(y))) {

  if(opt_level == 1) {
    lam = lam_opt(y, alpha, eps0, lam0)
    eps = eps_opt(y, alpha, lam, eps0)
  } else if (opt_level == 2) {
    theta = param_opt(y, alpha, eps0, lam0)
    lam = theta[1]
    eps = theta[2]
  } else {
    lam = lam0
    eps = eps0
  }

  yhat = tv1d(y, lam)
  yhat = l1fit(y, which(diff(yhat) != 0))

  # reduce cp
  cp_sets = cplink(yhat, eps)
  cp = cpfilt(y, yhat, cp_sets, alpha)

  return(cp)
}

#' @export
fmci = function(h, alpha = 1e-3, opt_level = 1, lam0 = sqrt(ncol(h)), eps0 = sqrt(ncol(h))) {
  arclen = ftv(h)
  arclen.cp = mci(arclen, alpha, opt_level, lam0, eps0)

  pc1 = fpc1(h)
  pc1.cp = mci(pc1, alpha, opt_level, lam0, eps0)

  round(poolC(sort(unique(c(arclen.cp, pc1.cp, ncol(h)))), sqrt(ncol(h))), 0)
}
