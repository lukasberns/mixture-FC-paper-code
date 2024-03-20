set.seed(1387517)

Ntrue = 16
theta.uniq.true = seq(-pi,pi,length=Ntrue+1)[-1]

Ntoys.per.true = 10000
Ntoys = Ntrue*Ntoys.per.true
theta.true = array(theta.uniq.true, Ntoys) # [Ntoys] but can reshape to [Ntrue,Ntoys.per.true]

# needs adjustment
epsilon.plot = 0.3

setwd(sprintf("20220528-mixturepaper-expfamily_rsrc/20240319-01-T2K-%dtrue-%dtoys", Ntrue, Ntoys.per.true))

get.lambda = function(theta) {
  B = 10
  Ntheta = length(theta)
  b = rep(seq(B), Ntheta)
  theta = rep(theta, each=B)
  A = 40.
  J = -1/4
  f = 1.
  epsilon = (b - (B+1)/2.) / B
  P = 0.25 * (1. - epsilon*epsilon)
  return(array(A * P * (1 + J * sin(theta + f*epsilon)), c(B,Ntheta)))
}

# observations
Nnu = length(get.lambda(0.))

lambda.true = get.lambda(theta.true)
nu.toy = array(rpois(Nnu*Ntoys, lambda.true), c(Nnu,Ntoys))

# now compute the chi2 at each POI value
Nfit = 100
theta.fit = seq(-pi,pi,length=Nfit+1)[-1]
lambda.fit = array(get.lambda(theta.fit), c(Nnu,Nfit))

get.chi2.pertoy = function(nu.toy, lambda.toy) {
  # @param  nu.toy [Nnu,Ntoys]
  # @param  lambda.toy [Nnu,Ntoys]
  # @return [Ntoys]
  
  nobs = nu.toy
  npred = lambda.toy
  
  nll = colSums((npred - nobs) - ifelse(nobs > 0., nobs*(log(npred) - log(nobs)), 0.))
  2.*nll
}

get.chi2 = function(nu.toy, lambda.fit) {
  # @param  nu.toy [Nnu,Ntoys]
  # @param  lambda.fit [Nnu,Nfit]
  # @return [Nfit,Ntoys]
  
  t(apply(lambda.fit, 2, function(lambda.one.fit) {
    # dpois(nu.toy, lambda.one.fit %o% rep(1., Ntoys), log=TRUE)
    nobs = nu.toy
    npred = lambda.one.fit %o% rep(1., dim(nu.toy)[2])
    get.chi2.pertoy(nobs, npred)
  }))
}

chi2.fit.toy = get.chi2(nu.toy, lambda.fit) # [Nfit, Ntoys]
# sum(is.na(chi2.fit.toy))

theta.toy = apply(chi2.fit.toy, 2, function(chi2.fit.onetoy) {
  theta.fit[which.min(chi2.fit.onetoy)]
}) # [Ntoys]
chi2.toy.min = apply(chi2.fit.toy, 2, function(chi2.fit.onetoy) {
  min(chi2.fit.onetoy)
}) # [Ntoys]

chi2.toy.true = get.chi2.pertoy(nu.toy, lambda.true) # [Ntoy]
Dchi2.toy.true = chi2.toy.true - chi2.toy.min # [Ntoys]

# now reformat
# sanity check:
assertthat::are_equal(array(theta.true,c(Ntrue,Ntoys.per.true))[,1], theta.uniq.true)
# then
Dchi2.toy.true = array(Dchi2.toy.true,c(Ntrue,Ntoys.per.true))




# critical values and binomial CI for conventional FC

Nprob = 5
qprobs = pchisq((1:Nprob)^2,1)

qprobs.ci = sapply(qprobs, function(p) { binom.test(round(p*Ntoys.per.true), Ntoys.per.true, conf.leve=0.683)$conf.int }) # [Nci,Nprobs]

Dchi2c.max = 100.
Dchi2c = apply(Dchi2.toy.true, 1, function(Dchi2.per.true) {
  quantile(c(Dchi2.per.true,Dchi2c.max), qprobs)
}) # [Nprob,Ntrue]
Dchi2c.ci = apply(Dchi2.toy.true, 1, function(Dchi2.per.true) {
  quantile(c(Dchi2.per.true,Dchi2c.max), qprobs.ci)
}) # [Nci*Nprob,Ntrue]
Dchi2c.ci = array(Dchi2c.ci, c(2,Nprob,Ntrue))





# mixture FC

lambda.uniq.true = get.lambda(theta.uniq.true)

# sampling probability from each true value
chi2smpl.toy.mixt.alltrue = get.chi2(nu.toy, lambda.uniq.true) # [Ntrue,Ntoys]
chi2smpl.toy.mixt = apply(chi2smpl.toy.mixt.alltrue, 2, function(chi2.fit.onetoy) {
  -2.*log(mean(exp(-0.5*chi2.fit.onetoy)))
}) # [Ntoys]

# shorthands
chi2smpl.fit.toy = chi2.fit.toy
Dchi2.fit.toy = chi2.fit.toy - rep(1,Nfit) %o% chi2.toy.min

# weights
weight.mixt = exp(-0.5*t(t(chi2smpl.fit.toy) - chi2smpl.toy.mixt)) # [Nfit,Ntoys]

# critical values
Dchi2c.mixt.v3 = sapply(1:Nfit, function(ifit.trgt) {
  ord = order(Dchi2.fit.toy[ifit.trgt,], decreasing=TRUE)
  y.srt = Dchi2.fit.toy[ifit.trgt,][ord]
  w.srt = weight.mixt[ifit.trgt,][ord]
  c.srt = cumsum(w.srt) / length(ord)
  # plot(c.srt, type='l', log='y', ylim=c(1e-6,1))
  # plot(c.srt, y.srt, type='l', log='x', xlim=c(1e-6,1), ylim=c(0.,25.))
  
  approx(c.srt, y.srt, 1.-qprobs)$y
}) # [Nprob,Nfit.trgt]

Dchi2c.mixt = Dchi2c.mixt.v3


# bootstrap error estimate
Nboot = 10
Dchi2c.mixt.boot = sapply(1:Nboot, function(iboot) {
  itoy.boot = sample(1:Ntoys, Ntoys, replace=TRUE)
  sapply(1:Nfit, function(ifit.trgt) {
    ord = order(Dchi2.fit.toy[ifit.trgt,itoy.boot], decreasing=TRUE)
    y.srt = Dchi2.fit.toy[ifit.trgt,itoy.boot][ord]
    w.srt = weight.mixt[ifit.trgt,itoy.boot][ord]
    c.srt = cumsum(w.srt) / length(ord)
    # plot(c.srt, type='l', log='y', ylim=c(1e-6,1))
    # plot(c.srt, y.srt, type='l', log='x', xlim=c(1e-6,1), ylim=c(0.,25.))

    approx(c.srt, y.srt, 1.-qprobs)$y
  }) # [Nprob,Nfit.trgt]
}) # [Nprob*Nfit.trgt,Nboot]
# 70 sec
Dchi2c.mixt.boot = array(Dchi2c.mixt.boot, c(Nprob,Nfit,Nboot))

Dchi2c.mixt.boot.mean = apply(Dchi2c.mixt.boot, c(1,2), mean) # [Nprob,Nfit]
Dchi2c.mixt.boot.sd   = apply(Dchi2c.mixt.boot, c(1,2), sd  )

Dchi2c.mixt.ci = rbind(
  Dchi2c.mixt.boot.mean + Dchi2c.mixt.boot.sd,
  Dchi2c.mixt.boot.mean - Dchi2c.mixt.boot.sd
)
Dchi2c.mixt.ci = array(Dchi2c.mixt.ci, c(Nprob,2,Nfit))
Dchi2c.mixt.ci = aperm(Dchi2c.mixt.ci, c(2,1,3)) # [Nci,Nprob,Nfit]



##############
# Make plots #
##############

customPDF = function(file) {
  pdf(file=file, 5, 5, family="serif")
  par(lwd=1.5, mar=c(5.1, 4.1+0.5, 4.1, 2.1))
}

weight.mixt.mean = apply(weight.mixt,1,mean)
weight.mixt.sd   = apply(weight.mixt,1,sd)
weight.mixt.mean.se = weight.mixt.sd / sqrt(apply(weight.mixt,1,length))

# sanity check (should be about 1)
customPDF("sum-of-weights.pdf")
plot(theta.fit, weight.mixt.mean, ylim=1 + 2.*c(-1,1)*max(weight.mixt.mean.se), type='l', lwd=2, xaxs='i', xlab=expression("Target "*theta[t]), ylab='Sample mean of weights')
polygon(c(theta.fit,rev(theta.fit)), c(
  weight.mixt.mean + weight.mixt.mean.se,
  rev(
    weight.mixt.mean - weight.mixt.mean.se
  )
),col=rgb(0,0,0,0.2))
abline(h=1, lty=2)
dev.off()


# comparison with binomial CI as bands
plot(theta.uniq.true, theta.uniq.true, 'n', ylim=c(0.,6.), ylab='sqrt[Dchi2c]')
for (iprob in 1:Nprob) {
  polygon(c(theta.uniq.true,rev(theta.uniq.true)), c(sqrt(Dchi2c.ci[1,iprob,]),rev(sqrt(Dchi2c.ci[2,iprob,]))), col=rgb(0,0,0,0.2))
  lines(theta.uniq.true, sqrt(Dchi2c[iprob,]), type='b')
}
for (iprob in 1:Nprob) {
  polygon(c(theta.fit,rev(theta.fit)), c(sqrt(Dchi2c.mixt.ci[1,iprob,]),rev(sqrt(Dchi2c.mixt.ci[2,iprob,]))), col=rgb(1,0,0,0.2))
  # lines(theta.fit, sqrt(Dchi2c.mixt[iprob,]), type='b', col=2, pch=4)
}




# error bands (mixture) vs. error bars (conventional)
ifit.match.true = sapply(theta.uniq.true, function(t) { which.min(abs(theta.fit - t)) })
for (interpolateMixt in c(TRUE, FALSE)) {
  customPDF(sprintf("critical-bands%s.pdf", ifelse(interpolateMixt, "", "-nointerp")))
  plot(theta.uniq.true, theta.uniq.true, 'n', ylim=c(0.,6.), ylab=expression(sqrt(Delta*chi[c]^2)), xlab=expression(theta), xlim=c(-pi,pi), xaxs='i', yaxs='i')
  for (iprob in 1:Nprob) {
    # polygon(c(theta.uniq.true,rev(theta.uniq.true)), c(sqrt(Dchi2c.ci[1,iprob,]),rev(sqrt(Dchi2c.ci[2,iprob,]))), col=rgb(0,0,0,0.2))
    xoffset = 0
    if (iprob == 4) { xoffset = -0.02; }
    if (iprob == 5) { xoffset = +0.02; }
    arrows(
      x0=theta.uniq.true+xoffset,
      y0=sqrt(Dchi2c.ci[1,iprob,]),
      x1=theta.uniq.true+xoffset,
      y1=sqrt(Dchi2c.ci[2,iprob,]),
      length=0.02,
      code=3,
      angle=90
    )
    points(theta.uniq.true+xoffset, sqrt(Dchi2c[iprob,]))
  }
  for (iprob in 1:Nprob) {
    if (interpolateMixt) {
      polygon(c(theta.fit,rev(theta.fit)), c(sqrt(Dchi2c.mixt.ci[1,iprob,]),rev(sqrt(Dchi2c.mixt.ci[2,iprob,]))), col=rgb(1,0,0), border='darkred')
      # lines(theta.fit, sqrt(Dchi2c.mixt[iprob,]), type='b', col=2, pch=4)
    }
    else {
      arrows(
        x0=theta.uniq.true,
        y0=sqrt(Dchi2c.mixt.ci[1,iprob,ifit.match.true]),
        x1=theta.uniq.true,
        y1=sqrt(Dchi2c.mixt.ci[2,iprob,ifit.match.true]),
        length=0.02,
        code=3,
        angle=90,
        col=2,
        lwd=1.5
      )
    }
  }
  abline(h=sqrt(qchisq(qprobs,1)), lty=2, col=rgb(0,0,0,0.3))
  # legend('bottomleft', c('Conventional FC', 'Mixture model', 'Wilks'), lty=c(1,1,2), col=c("black","red",rgb(0,0,0,0.3)),pch=c(1,NA,NA), lwd=c(1,1.5,1), bty='n', inset=0.03, ncol=3)
  legend('bottom', c('Conventional FC', 'Mixture model'), lty=c(1,1), col=c("black","red"),pch=c(1,NA,NA), lwd=c(2,2.5), bty='n', ncol=2, x.intersp=0.6, text.width=2.2)
  dev.off()
}

# error bands (mixture) vs. error bars (conventional)
plot.crit.ylim.min = c(0.8, 1.8, 2.7, 3.1, 4.2)
plot.crit.ylim.max = c(1.3, 2.3, 3.2, 4.2, 5.2)
for (interpolateMixt in c(TRUE,FALSE)) {
  for (isig in seq(length(plot.crit.ylim.min))) {
    # zoom of 1 sigma critical values
    customPDF(sprintf("critical-bands%s-%dsig.pdf", ifelse(interpolateMixt, "", "-nointerp"), isig))
    plot(theta.uniq.true, theta.uniq.true, 'n', ylim=c(plot.crit.ylim.min[isig],plot.crit.ylim.max[isig]), ylab=expression(sqrt(Delta*chi[c]^2)), xlab=expression(theta), xlim=c(-pi,pi), xaxs='i', yaxs='i')
    abline(h=sqrt(qchisq(qprobs,1)), lty=2, col=rgb(0,0,0,0.3))
    for (iprob in isig) {
      polygon(
        c(
          c(theta.uniq.true-2*pi, theta.uniq.true),
          rev(c(theta.uniq.true-2*pi, theta.uniq.true))
        ),
        c(
          sqrt(rep(Dchi2c.ci[1,iprob,],2)),
          rev(sqrt(rep(Dchi2c.ci[2,iprob,],2)))
        ), col=rgb(0,0,0,0.15), border=FALSE)
   }
    for (iprob in isig) {
      if (interpolateMixt) {
        polygon(c(theta.fit,rev(theta.fit)), c(sqrt(Dchi2c.mixt.ci[1,iprob,]),rev(sqrt(Dchi2c.mixt.ci[2,iprob,]))), col=rgb(1,0,0), border='darkred')
        # lines(theta.fit, sqrt(Dchi2c.mixt[iprob,]), type='b', col=2, pch=4)
      }
      else {
        arrows(
          x0=theta.uniq.true,
          y0=sqrt(Dchi2c.mixt.ci[1,iprob,ifit.match.true]),
          x1=theta.uniq.true,
          y1=sqrt(Dchi2c.mixt.ci[2,iprob,ifit.match.true]),
          length=0.02,
          code=3,
          angle=90,
          col=2,
          lwd=1.5
        )
      }
    }
    for (iprob in isig) {
      # polygon(c(theta.uniq.true,rev(theta.uniq.true)), c(sqrt(Dchi2c.ci[1,iprob,]),rev(sqrt(Dchi2c.ci[2,iprob,]))), col=rgb(0,0,0,0.2))
      xoffset = 0
      if (iprob == 4) { xoffset = -0.02; }
      if (iprob == 5) { xoffset = +0.02; }
      arrows(
        x0=theta.uniq.true+xoffset,
        y0=sqrt(Dchi2c.ci[1,iprob,]),
        x1=theta.uniq.true+xoffset,
        y1=sqrt(Dchi2c.ci[2,iprob,]),
        length=0.02,
        code=3,
        angle=90
      )
      points(theta.uniq.true+xoffset, sqrt(Dchi2c[iprob,]))
    }
    # legend('bottomleft', c('Conventional FC', 'Mixture model', 'Wilks'), lty=c(1,1,2), col=c("black","red",rgb(0,0,0,0.3)),pch=c(1,NA,NA), lwd=c(1,1.5,1), bty='n', inset=0.03, ncol=3)
    legend('bottom', c('Conventional FC', 'Mixture model'), lty=c(1,1), col=c("black","red"),pch=c(1,NA,NA), lwd=c(2,2.5), bty='n', ncol=2, x.intersp=0.6, text.width=2.2)
    dev.off()
  }
}


# toy distribution
require(plotrix) # for weighted histogram

Dchi2.hist.bins = seq(0.,32.,len=100)
Dchi2.hist.bins.centers = (Dchi2.hist.bins[-1]+Dchi2.hist.bins[-length(Dchi2.hist.bins)])/2.

# the following two are -pi/2
itrue = 4
ifit.trgt = 25
theta.uniq.true[itrue]
theta.fit[ifit.trgt]
C.plot = 1.

# more toys as ground truth for conventional method
# this block takes a few min to compute (see gt.Ntoys.per.true below)
gt.Ntrue = 1
gt.theta.uniq.true = theta.uniq.true[itrue]
gt.Ntoys.per.true = 10000000 # 1M = 2 min, 800 MB. 4M = 8 min, 3.2 GB. 10M
gt.Ntoys = gt.Ntrue*gt.Ntoys.per.true
gt.theta.true = array(gt.theta.uniq.true, gt.Ntoys) # [Ntoys] but can reshape to [Ntrue,Ntoys.per.true]
# observations
gt.lambda.true = get.lambda(gt.theta.true)
gt.nu.toy = array(rpois(Nnu*gt.Ntoys, gt.lambda.true), c(Nnu,gt.Ntoys))

# since we will be limited by memory, we will process in batches
gt.Ntoy.batches = c(seq(1, gt.Ntoys.per.true, 1e5), gt.Ntoys.per.true+1)
gt.chi2.toy.min = sapply(seq(gt.Ntoy.batches[-1]), function(ibatch) {
  print(sprintf("%3d / %3d", ibatch, length(gt.Ntoy.batches)-1))
  gt.itoy = seq(gt.Ntoy.batches[ibatch], gt.Ntoy.batches[ibatch+1]-1)
  gt.chi2.fit.toy = get.chi2(gt.nu.toy[,gt.itoy], lambda.fit) # [Nfit, gt.Ntoys]
  gt.chi2.toy.min.here = apply(gt.chi2.fit.toy, 2, function(chi2.fit.onetoy) {
    min(chi2.fit.onetoy)
  }) # [gt.Ntoys]
  remove(gt.chi2.fit.toy) # this object is 800 MB for 1M toys
  return(gt.chi2.toy.min.here)
}) # [Ntoys per batch, Nbatch]
gt.chi2.toy.min = c(gt.chi2.toy.min) # squash, this orders them as c(Ntoys per batch1, Ntoys per batch2, ...)
# gt.chi2.fit.toy = get.chi2(gt.nu.toy, lambda.fit) # [Nfit, gt.Ntoys]
# # sum(is.na(gt.chi2.fit.toy))
# gt.chi2.toy.min = apply(gt.chi2.fit.toy, 2, function(chi2.fit.onetoy) {
#   min(chi2.fit.onetoy)
# }) # [gt.Ntoys]
# remove(gt.chi2.fit.toy) # this object is 800 MB

gt.chi2.toy.true = get.chi2.pertoy(gt.nu.toy, gt.lambda.true) # [gt.Ntoy]
gt.Dchi2.toy.true = gt.chi2.toy.true - gt.chi2.toy.min # [gt.Ntoys]
# now reformat
# sanity check:
assertthat::are_equal(array(gt.theta.true,c(gt.Ntrue,gt.Ntoys.per.true))[,1], gt.theta.uniq.true)
# then
gt.Dchi2.toy.true = array(gt.Dchi2.toy.true,c(gt.Ntrue,gt.Ntoys.per.true))


Dchi2.hist.std = weighted.hist(
  pmax(pmin(c(Dchi2.toy.true[itrue,]),max(Dchi2.hist.bins)), 0.),
  breaks=Dchi2.hist.bins,
  plot=FALSE
)$counts
Dchi2.hist.std.nonneg = ifelse(Dchi2.hist.std > 0, Dchi2.hist.std, 1e-11)

gt.Dchi2.hist.std = weighted.hist(
  pmax(pmin(c(gt.Dchi2.toy.true[1,]),max(Dchi2.hist.bins)), 0.),
  breaks=Dchi2.hist.bins,
  plot=FALSE
)$counts
gt.Dchi2.hist.std.nonneg = ifelse(gt.Dchi2.hist.std > 0, gt.Dchi2.hist.std, 1e-11)

Dchi2.hist.mixt.unweighted = weighted.hist(
  pmax(pmin(c(Dchi2.fit.toy[ifit.trgt,]),max(Dchi2.hist.bins)), 0.),
  breaks=Dchi2.hist.bins,
  plot=FALSE
)$counts
Dchi2.hist.mixt.unweighted.nonneg = ifelse(Dchi2.hist.mixt.unweighted > 0, Dchi2.hist.mixt.unweighted, 1e-11)

Dchi2.hist.mixt = weighted.hist(
  pmax(pmin(c(Dchi2.fit.toy[ifit.trgt,]),max(Dchi2.hist.bins)), 0.),
  weight.mixt[ifit.trgt,],
  breaks=Dchi2.hist.bins,
  plot=FALSE
)$counts/Ntrue
Dchi2.hist.mixt.w2 = weighted.hist(
  pmax(pmin(c(Dchi2.fit.toy[ifit.trgt,]),max(Dchi2.hist.bins)), 0.),
  weight.mixt[ifit.trgt,]^2,
  breaks=Dchi2.hist.bins,
  plot=FALSE
)$counts/Ntrue^2
Dchi2.hist.mixt.nonneg    = ifelse(Dchi2.hist.mixt > 0, Dchi2.hist.mixt, 1e-11)
Dchi2.hist.mixt.w2.nonneg = ifelse(Dchi2.hist.mixt > 0, Dchi2.hist.mixt.w2, 1e-11)
Dchi2.hist.mixt.se = sqrt(Dchi2.hist.mixt.w2.nonneg)

histWithoutErrorBars = function(bins, n, col=1, lty=1) {
  # bins will have one more entry than n
  lines(bins, c(n,1e-11), type='s', col=col, lty=lty)
}
histWithBinomErrorBars = function(binCenters, n, col=1, showZeros=FALSE, pch=20, scale=1., lwd=2, lend="round") {
  # ci = binconf(n, sum(n), alpha=1-pchisq(1^2,1), "wilson")
  # ci = BinomCI(n, sum(n), conf.level = pchisq(1^2,1), method = "wilson")
  ci = sapply(n, function(x) { binom.test(x, sum(n), conf.level=pchisq(1^2,1))$conf.int*sum(n) }) # clopper-pearson
  arrows(
    x0=binCenters[showZeros | n>0],
    y0=pmax(scale*ci[1,],1e-11)[showZeros | n>0],
    x1=binCenters[showZeros | n>0],
    y1=(scale*ci[2,])[showZeros | n>0],
    length=0.,
    code=3,
    angle=90,
    col=col,
    lwd=lwd,
    lend=lend
  )
  hollowCircle = 1
  solidCircle = 20
  points(binCenters, n*scale, col=col, pch=pch)
}
histWithErrorBands = function(bins, w, w2, col=1, border=1) {
  # bins will have one more entry than w and w2
  se = sqrt(w2)
  polygon(
    c(rep(bins,each=2)[-1],rev(rep(bins,each=2)[-1])),
    c(rep(w+se,each=2),1e-11,1e-11,rev(rep(w-se,each=2))),
    col=col,
    border=border
  )
}

# plot
for (withMixture in c(FALSE, TRUE)) {
  customPDF(sprintf("toys-std%d-mixture%d%s.pdf", itrue, ifit.trgt, ifelse(withMixture, "", "-nomix")))
  plot( Dchi2.hist.bins, pmax(Dchi2.hist.bins,1e-11), type='n', col=2, xaxs='i', log='y', ylim=range(subset(data.frame(v=c(Dchi2.hist.mixt,Dchi2.hist.std,Dchi2.hist.mixt.unweighted)), v>0)$v), ylab='Number of pseudo-experiments', xlab=expression(Delta*chi[t]^2), xlim=c(0,25))
  
  for (iprob in 1:Nprob) {
    abline(v=Dchi2c.mixt[iprob,ifit.trgt], col='gray', lty=2)
  }
  
  histWithBinomErrorBars(Dchi2.hist.bins.centers, Dchi2.hist.std, col=1, showZeros=TRUE)
  histWithBinomErrorBars(Dchi2.hist.bins.centers, gt.Dchi2.hist.std, col=1, showZeros=TRUE, lwd=4, lend="butt", pch=NA, scale=Ntoys.per.true/gt.Ntoys)
  histWithBinomErrorBars(Dchi2.hist.bins.centers, Dchi2.hist.mixt.unweighted, col=4, showZeros=TRUE, pch=18)
  
  # histWithoutErrorBars(Dchi2.hist.bins, Dchi2.hist.mixt, col=2)
  if (withMixture) {
    histWithErrorBands(Dchi2.hist.bins, Dchi2.hist.mixt.nonneg, Dchi2.hist.mixt.w2.nonneg, col=rgb(1,0,0,0.5), border=2)
  }
  # lines(Dchi2.hist.bins.centers, Dchi2.hist.mixt+Dchi2.hist.mixt.se, type='s', col=5)
  # lines(Dchi2.hist.bins.centers, Dchi2.hist.mixt-Dchi2.hist.mixt.se, type='s', col=5)
  
  histWithoutErrorBars(Dchi2.hist.bins, Dchi2.hist.mixt.nonneg * exp(+0.5*Dchi2.hist.bins.centers-0.5*epsilon.plot), col=3, lty=2)
  # histWithoutErrorBars(Dchi2.hist.bins, Dchi2.hist.std.nonneg, col=1)
  # histWithoutErrorBars(Dchi2.hist.bins, Dchi2.hist.mixt.unweighted.nonneg, col=4)
  
  leg.items = c('Conventional FC',sprintf('Conv. FC (%.0fx, scaled)', gt.Ntoys/Ntoys.per.true),'Weighted mixture','Unweighted mixture','Theoretical lower limit')
  leg.col = c(1,1,2,4,3)
  leg.pch = c(20,NA,NA,18,NA)
  leg.lty = c(1,1,1,1,2)
  leg.lwd = c(2,4,2,2,2)
  if (!withMixture) {
    leg.items = leg.items[-3]
    leg.col = leg.col[-3]
    leg.pch = leg.pch[-3]
    leg.lty = leg.lty[-3]
    leg.lwd = leg.lwd[-3]
  }
  # legend('bottomleft', leg.items,lty=leg.lty,pch=leg.pch,col=leg.col,bg="white")
  legend('bottomleft', leg.items,lty=leg.lty,col=leg.col,lwd=leg.lwd,bg="white")
  legend('bottomleft', rep("",length(leg.items)),pch=leg.pch,col=leg.col,bty='n',inset=c(0.04,0.))
  dev.off()
}

# scaled so we can see consistency within the error bars
# we choose as reference the Gaussian best-fit among the 1000x convFC (better at low Dchi2) and mixture result (better at high Dchi2)
ref.hist.a1 = Dchi2.hist.mixt.nonneg
ref.hist.a2 = gt.Dchi2.hist.std * (Ntoys.per.true/gt.Ntoys)
ref.hist.w1 = 1./Dchi2.hist.mixt.w2.nonneg
ref.hist.w2 = 1./(gt.Dchi2.hist.std+1.) / (Ntoys.per.true/gt.Ntoys)**2
ref.hist = (ref.hist.a1*ref.hist.w1 + ref.hist.a2*ref.hist.w2) / (ref.hist.w1 + ref.hist.w2)

plot( Dchi2.hist.bins, pmax(Dchi2.hist.bins,1e-11), type='n', col=2, xaxs='i', ylim=c(0.9, 1.1), ylab='Number of pseudo-experiments / reference', xlab=expression(Delta*chi[t]^2), xlim=c(0,25))
for (iprob in 1:Nprob) {
  abline(v=Dchi2c.mixt[iprob,ifit.trgt], col='gray', lty=2)
}
# histWithBinomErrorBars(Dchi2.hist.bins.centers, Dchi2.hist.std, col=1, showZeros=TRUE, scale=1./ref.hist)
histWithBinomErrorBars(Dchi2.hist.bins.centers, gt.Dchi2.hist.std, col=1, showZeros=TRUE, lwd=4, lend="butt", pch=NA, scale=Ntoys.per.true/gt.Ntoys/ref.hist)
histWithErrorBands(Dchi2.hist.bins, Dchi2.hist.mixt.nonneg/ref.hist, Dchi2.hist.mixt.w2.nonneg/ref.hist/ref.hist, col=rgb(1,0,0,0.5), border=2)


# chi2 for agreement.
# difference between the two estimates
# and normalize by the error on this difference
# assuming Gaussianity.

ref.s = gt.Ntoys/Ntoys.per.true
ref.chi2consistency = (ref.hist.a1 - ref.hist.a2)**2 / (Dchi2.hist.mixt.w2.nonneg + gt.Dchi2.hist.std/ref.s**2)
# cut at 20 so we are not affected by zero toys for conv FC
sum(ref.chi2consistency[Dchi2.hist.bins <= 20.])
# 65.12502
sum(Dchi2.hist.bins <= 20.)
# 62
1. - pchisq(sum(ref.chi2consistency[Dchi2.hist.bins <= 20.]), sum(Dchi2.hist.bins <= 20.))
# 0.3685405


# plot of estimated errors on critical values
Pmixt = Dchi2.hist.mixt.nonneg / Ntoys.per.true
Dchi2max.plot = rev(Dchi2.hist.bins)[1]
Pmax = rev(Pmixt)[1] # technically not exact but should do

y.plot = Dchi2.hist.bins.centers
ymax.plot = pmax(y.plot,Dchi2max.plot)
A.plot = 1/exp(0.5*(y.plot - epsilon.plot))
B.plot = 1. - A.plot
S.plot = Ntrue
gamma.plot = (A.plot + B.plot*Pmax/Pmixt - 1/S.plot*Pmixt) / (1 - Pmixt)

dPconv = sqrt((1.-Pmixt)/Pmixt/sum(Dchi2.hist.std))
dPmixt.upper = dPconv * sqrt(gamma.plot)

customPDF(sprintf("relerr-std%d-mixture%d.pdf", itrue, ifit.trgt))
plot(Dchi2.hist.bins.centers, Dchi2.hist.bins.centers, 'n', col=1, ylim=c(1e-2,0.2), log='y', ylab="Estimated relative error on CDF estimate", xlab=expression(Delta*chi[t]^2), xaxs='i', xlim=c(0,25))
lines(Dchi2.hist.bins.centers, dPconv, 'l', col=1, lty=5)
lines(Dchi2.hist.bins.centers, dPmixt.upper, 'l', col=4, lty=3)
lines(Dchi2.hist.bins.centers, Dchi2.hist.mixt.se / Dchi2.hist.mixt.nonneg, 'l', col=2)
legend('bottomright', c('Conventional FC','Mixture FC (bootstrap)',expression('Upper bound on '*gamma)),lty=c(5,1,3),col=c(1,2,4),bty="n")
dev.off()




# toys at various theta
for (withMixture in c(FALSE, TRUE)) {
  customPDF(sprintf("toys-fewtrue%s.pdf", ifelse(withMixture, "", "-nomix")))
  plot( Dchi2.hist.bins, 0.*Dchi2.hist.bins, type='n', col=2, xaxs='i', log='y', yaxs='i', ylim=c(0.8,1e5), ylab='Number of pseudo-experiments', xlab=expression(Delta*chi[t]^2), xlim=c(0,25))

  for (iprob in 1:Nprob) {
    abline(v=Dchi2c.mixt[iprob,ifit.trgt], col='gray', lty=2)
  }

  df.pts.name =  c(expression("Sampling "*theta*" = "*-pi*"/2 (target)"), expression("Sampling "*theta*" = 0"), expression("Sampling "*theta*" = "*+pi*"/2"), bquote("Mixture of "*.(sprintf("%s", Ntrue))*" "*theta*"-values"))
  df.pts = data.frame(
    col = c(1,2,3,4),
    pch = c(20,15,17,18),
    itrue = c( Ntrue%/%4, Ntrue %/% 2, 3*Ntrue%/%4, -1)
  )
  
  if (!withMixture) {
    df.pts = df.pts[-length(df.pts.name),]
    df.pts.name = df.pts.name[-length(df.pts.name)]
  }

  for (ipt in seq(df.pts$col)) {
    itrue = df.pts$itrue[ipt]
    ifit.trgt = 25 # always -pi/2

    d = Dchi2.fit.toy[ifit.trgt,]
    if (itrue > 0) {
      d = array(d, c(Ntrue,Ntoys.per.true))[itrue,]
    }
    Dchi2.hist.mixt.unweighted = weighted.hist(
      pmax(pmin(c(d),max(Dchi2.hist.bins)), 0.),
      breaks=Dchi2.hist.bins,
      plot=FALSE
    )$counts
    Dchi2.hist.mixt.unweighted.nonneg = ifelse(Dchi2.hist.mixt.unweighted > 0, Dchi2.hist.mixt.unweighted, 1e-11)

    # lines(Dchi2.hist.bins.centers, Dchi2.hist.mixt.unweighted.nonneg , type='s', col=df.pts$col[ipt])
    # histWithoutErrorBars(Dchi2.hist.bins, Dchi2.hist.mixt.unweighted.nonneg, col=df.pts$col[ipt])
    histWithBinomErrorBars(Dchi2.hist.bins.centers, Dchi2.hist.mixt.unweighted, col=df.pts$col[ipt], pch=df.pts$pch[ipt])
  }
  legend('topright', legend=df.pts.name, lty=1, col=df.pts$col, pch=df.pts$pch, box.col="white", inset=c(0.005,0.005))
  dev.off()
}

# epsilon
chi2.toy.mingrid = apply(chi2smpl.toy.mixt.alltrue, 2, function(chi2.fit.onetoy) {
  min(chi2.fit.onetoy)
}) # [Ntoys]

# the following two are -pi/2
ifit.trgt = 25
theta.fit[ifit.trgt]

customPDF(sprintf("epsilon.pdf", itrue, ifit.trgt))
plot(
  Dchi2.fit.toy[ifit.trgt,1:10000],
  (chi2.toy.mingrid - chi2.toy.min)[1:10000],
  pch='.',
  xaxs='i',
  yaxs='i',
  ylim=c(0.,1.1*epsilon.plot),
  col=2,
  xlab=expression(Delta*chi[t]^2),
  ylab=expression(Delta*chi^2*(hat(theta)[S]))
)
abline(h=epsilon.plot, lty=2)
dev.off()




# the following two are -pi/2
itrue = 4
ifit.trgt = 25
theta.uniq.true[itrue]
theta.fit[ifit.trgt]

Dchi2.plot = seq(0., 40., len=400)
wL.plot = exp(-0.5*Dchi2.plot)
wU.plot = Nnu * exp(-0.5*(Dchi2.plot - epsilon.plot))

customPDF(sprintf("weights.pdf", itrue, ifit.trgt))
plot(
  Dchi2.fit.toy[ifit.trgt,1:10000],
  weight.mixt[ifit.trgt,1:10000],
  pch='.',
  log='y',
  xaxs='i',
  col=2,
  xlab=expression(Delta*chi[t]^2),
  ylab='Pseudo-experiment weight',
  xlim=range(Dchi2.plot)
)
abline(h=1, lty=4)
lines(Dchi2.plot, wL.plot, lty=2, col=4)
lines(Dchi2.plot, wU.plot, lty=2, col=4)
legend('bottomleft', c('Conventional Feldman-Cousins','Mixture Feldman-Cousins','Theoretical bounds'), col=c(1,2,4), lty=c(4,NA,2), pch=c(NA,20,NA), pt.cex=0.5, bty='n', inset=c(0.05,0))
dev.off()



# Plot upper bound on gamma

setwd(sprintf("../20240319-02-generic"))

customPDF = function(file) {
  pdf(file=file, 5, 5, family="serif")
  par(lwd=1.5, mar=c(5.1, 4.1+0.5, 4.1, 2.1))
}

customPDF("gamma.pdf")
Dchi2max = 35
y = seq(0,40,len=400)
ymax = pmax(y,Dchi2max)
k = 1
P = 1.-pchisq(y,k)
Pmax = 1.-pchisq(ymax,k)

epsilon = 1
A = 1/exp(0.5*(y - epsilon))
B = 1 - A
S = 10
gamma = (A + B*Pmax/P - 1/S*P) / (1 - P)
plot(y, A/(1-P), 'l', log='y', col=2, lty=2, ylab=expression('Upper bound on '*gamma),
ylim=c(0.5e-4,3*S), xaxs='i', yaxs='i', yaxt="n") # , xaxt="n")
lines(y, B*Pmax/P/(1-P), col=3, lty=4)
lines(y, gamma, col=1)
abline(h=1, lty=3)
axis(2, at=c(1,1e-2,1e-4), labels=c(expression(1), expression(10^-2), expression(10^-4)), las=1)
axis(1, at=c(epsilon, Dchi2max), labels=c(expression(epsilon), expression(Delta*chi[max]^2)), col.axis=4, col=4)
abline(v=c(epsilon,Dchi2max), lty=2, col=4)
legend("top", c("Total", expression(y <= Y(x)*" < "*Delta*chi[max]^2), expression(Y(x) >= Delta*chi[max]^2)), lty=c(1,2,4), col=1:3, bty="n")
box()
dev.off()


# Plot relative error on estimated CDF

customPDF("variance.pdf")
nexp = 10000
varCDFconv = 1/nexp*P*(1-P)
varCDFmix  = varCDFconv*gamma

plot(y, sqrt(varCDFconv)/P, 'l', log='y', col=1, lty=4, ylab=expression('Relative error on estimated CDF'), xaxs='i', yaxs='i', yaxt="n")
lines(y, sqrt(varCDFmix)/P, col=2)
axis(2, at=c(1e-2,1e-1,1,1e1,1e2), labels=c(expression(10^-2), expression(10^-1), 1, 10, 100), las=1)
axis(1, at=c(epsilon, Dchi2max), labels=c(expression(epsilon), expression(Delta*chi[max]^2)), col.axis=4, col=4)
abline(h=0.1, lty=3)
abline(v=c(epsilon,Dchi2max), lty=2, col=4)
legend("top", c("Conventional FC", "Mixture FC (upper bound)"), lty=c(4,1), col=1:2, bty="n")
box()
dev.off()



# Plot expected number of events for different values of dCP #

customPDF("lambda.pdf")
plot(seq(Nnu+1), type='n', xaxs='i', yaxs='i', ylim=c(0,14), xlab="Bin index", ylab="Predicted events per bin")
histWithoutErrorBars(seq(0,Nnu)+1, get.lambda(-pi/2), col=3, lty=6)
histWithoutErrorBars(seq(0,Nnu)+1, get.lambda(0)    , col=1, lty=1)
histWithoutErrorBars(seq(0,Nnu)+1, get.lambda(+pi/2), col=4, lty=4)
histWithoutErrorBars(seq(0,Nnu)+1, get.lambda(pi)   , col=2, lty=5)
legend("bottomleft", c(expression(delta[CP] == -pi/2),
                   expression(delta[CP] ==  0),
                   expression(delta[CP] == +pi/2),
                   expression(delta[CP] ==  pi)), col=c(3,1,4,2), lty=c(6,1,4,5), bty="n")
box()
dev.off()