#' Simulation for functional quadratic regression

# generate x
library(fda)
library(MASS)
source("Functions_freg.r")

m_ni = 21   # number of time points
m_ti = seq(0, 1, length=m_ni)  #time points (equally spaced)

# simulation setting for true functions x
nbasis = 7
bbasis0 = create.bspline.basis(c(0,1), nbasis = nbasis) 
bbasis02 = eval.basis(m_ti, bbasis0)  #basis matrix: m_ti x nbasis
JPhi0 = inprod(bbasis0, bbasis0)

# coefficients for beta
sig_beta = toeplitz(round((nbasis:1)/nbasis, 3)) # variance of Wishart dist.
coeff_beta = rWishart(1, 10, sig_beta)   # generate coefficent of coefficient function by Wishart dist.
beta = bbasis02 %*% coeff_beta[,,1] %*% t(bbasis02)  #m_ti x m_ti

# coefficients for gamma
sig_gamma = toeplitz(round((nbasis:1)/nbasis, 3)) # variance of Wishart dist.
coeff_gamma = rWishart(nbasis, 10, sig_gamma) # generate coefficent of coefficient function by Wishart dist.
tbuf = matrix(nr=nbasis^2, nc=nbasis)
for(i in 1:nbasis)
  tbuf[(1:nbasis)+(i-1)*nbasis, ] = coeff_gamma[,,i]
gamma = (bbasis02 %x% bbasis02) %*% (tbuf) %*% t(bbasis02)  #m_ti^2 x m_ti éQçlÇ∆ÇµÇƒåvéZ




#########################################################################################
library(fda)
library(MASS)
set.seed(Sys.time())
m_n = 50   # sample size
m_t = rep(1:m_ni, m_n)  # time point number for all subjects
m_t1 = c(1, (1:m_n)*m_ni+1)  #start point of m_t for each subject
m_sig_GP = 0.6  #noise level for GP
m_sig_errx = 0.3  #noise level for x
m_sig_erry = 0.3  #noise level for y

# random coefficient for x
mu_x = rep(0, m_n)
sig_x = toeplitz(round(0.9^(-0:(m_n-1)), 3))

# response function
yt = matrix(nr=m_n, nc=m_ni)  #true y
yf = yt  # y + noise process
y = yt   # y + noise process + noise

# for estimation in iteration
m_mx = 6
m_my = 6
bbasis = create.bspline.basis(c(0,1), nbasis=m_mx)
bbasis2 = eval.basis(m_ti, bbasis)  #basis matrix: m_ti x nbasis
JPhi = inprod(bbasis, bbasis)
quadcoeff = matrix(nr=m_mx^2, nc=m_n)  #coefficients for quadratic term


########################################################################################
#iteration start
iter.max = 100
# result
icnum = 4
resMSE = matrix(nr=iter.max, nc=icnum+6)
resLambda = matrix(nr=iter.max, nc=icnum)
for(iter in 1:iter.max){
    
  # generate random process x
  coeff_x = mvrnorm(nbasis, mu_x, 0.5*sig_x)
  xf = t(coeff_x) %*% t(bbasis02)
  x = xf  # xf + noise
  
  # check random process x
  if(0==1){
    plot(xf[1,],type="l", ylim=range(xf))
    for(i in 2:m_n) lines(xf[i,])
  }

  # generate functions y and yf
  for(i in 1:m_n){
    yt[i, ] = coeff_x[, i] %*% JPhi0 %*% coeff_beta[,,1] %*% t(bbasis02) +  #linear order
              (coeff_x[, i] %x% coeff_x[,i]) %*% (JPhi0 %x% JPhi0) %*% tbuf %*% t(bbasis02) #quadratic order
    yf[i, ] = yt[i,] + GP(m_ni, r=m_sig_GP)
  }
  
  # check functions y
  if(0==1){
    plot(yt[1,],type="l", ylim=range(yt))
    for(i in 2:m_n) lines(yt[i,])
  }
  # add noise  to  x and y
  for(i in 1:m_n){
    rngy = diff(range(yf[i,]))
    y[i,] = yf[i,] + rnorm(m_ti, 0, m_sig_erry*rngy)
    rngx = diff(range(xf[i,]))
    x[i,] = xf[i,] + rnorm(m_ti, 0, m_sig_errx*rngx)
  }

  # check functions y
  if(0==1){
    plot(y[1,],type="l", ylim=range(y))
    for(i in 2:m_n) lines(y[i,])
  }
  
    
  ###### data generation end ######

  ## transform x and y into functional data
  argvals = matrix(m_ti, nr=m_n, nc=m_ni, byrow=T)
  wmat = smooth.basis(argvals=t(argvals), y=t(x), fdParobj=bbasis)$fd$coefs

  
  #estimate and calculate prediction error
  # Freg-quad-interaction
  # create design matrix
  for(i in 1:m_n)
    quadcoeff[, i] = (JPhi %x% JPhi) %*% (wmat[, i] %x% wmat[, i])
  Z = t(rbind(rep(1, m_n), JPhi %*% wmat, quadcoeff)) # ncol=1 + m_x + m_x^2
  
  # Penalized likelihood method
  minic = rep(Inf, icnum) 
  for(i in 1:8){
    lambda = 10^(-i/1.5-1)
    est = try(CalcPGPR(Z, y, bbasis2, lambda))
    if(!is.null(attr(est,"class"))) next
    ic = IC_PGPR(Z, y, bbasis2, est, lambda)
    for(j in 1:icnum){
      if(ic[j] < minic[j]){
        minic[j] = ic[j]
        resLambda[iter, j] = lambda
      }
    }
  }
  for(j in 1:icnum){
    est = try(CalcPGPR(Z, y, bbasis2, resLambda[iter,j]))
    yhat = Z %*% t(est$Theta) %*% t(bbasis2)
    resMSE[iter, j+6] = sum((yt-yhat)^2)
  }
  
  #MLE
  est = try(CalcGPR(Z, y, bbasis2))
  if(is.null(attr(est,"class"))){
    yhat = Z %*% t(est$Theta) %*% t(bbasis2)
    resMSE[iter, 1] = sum((yt-yhat)^2)
  }
  
  # SFreg-interaction
  est = try(ginv(t(Z)%*%Z)%*%t(Z)%*%y)
  if(is.null(attr(est,"class"))){
    yhat = Z%*%est
    resMSE[iter, 3] = sum((yt-yhat)^2)
  }
  ######################

  # Freg-linear
  Z = t(rbind(rep(1, m_n), JPhi %*% wmat))
  est = try(CalcGPR(Z, y, bbasis2))
  if(is.null(attr(est,"class"))){
    yhat = Z %*% t(est$Theta) %*% t(bbasis2)
    resMSE[iter, 2] = sum((yt-yhat)^2)
  }
  
  
  # Sreg-quad
  Z = t(rbind(rep(1, m_n), t(x), t(x)^2))
  est = try(ginv(t(Z)%*%Z)%*%t(Z)%*%y)
  if(is.null(attr(est,"class"))){
    yhat = Z%*%est
    resMSE[iter, 4] = sum((yt-yhat)^2)
  }
  # Sreg-quad-interaction
  Z = t(rbind(rep(1, m_n), t(x)))
  for(j in 1:m_ni){
    for(k in 1:j){
      Z = cbind(Z, x[, j]*x[, k])
    }
  }
  est = try(ginv(t(Z)%*%Z)%*%t(Z)%*%y)
  if(is.null(attr(est,"class"))){
    yhat = Z%*%est
    resMSE[iter, 5] = sum((yt-yhat)^2)
  }
  # Sreg-linear
  Z = t(rbind(rep(1, m_n), t(x)))
  est = try(ginv(t(Z)%*%Z)%*%t(Z)%*%y)
  if(is.null(attr(est,"class"))){
    yhat = Z%*%est
    resMSE[iter, 6] = sum((yt-yhat)^2)  
  }

  print(iter)
}

#resMSE
apply(resMSE / (m_n*m_ni),2,mean)
apply(resMSE / (m_n*m_ni),2,sd)

boxplot(resMSE/(m_n*m_ni))

#plot predicted values
plot(yt[1,],ylim=range(yt),type="n")
for(i in 1:m_n){
#  plot(yt[1,],ylim=range(yt),type="n")
  lines(yt[i,],col=i)
  lines(yhat[i,],col=i,lty=2)
#  Sys.sleep(0.1)
}



#boxplot
# MLE
library(ggplot2)
x <- data.frame(
  x = as.factor(c(rep("FF-INT",iter.max),rep("FF-LIN",iter.max),rep("SF-INT",iter.max),
                  rep("INT",iter.max),rep("QUAD",iter.max),rep("LIN",iter.max))),
  y = c(resMSE[,1],resMSE[,2],resMSE[,3],resMSE[,5],resMSE[,4],resMSE[,6]) / (m_n*m_ni),
  z=1:600
)
g <- ggplot(x,aes(x=reorder(x,z),y=y), xlab="")
g <- g + geom_boxplot() + xlab("") + ylab("ASE") +
  theme_bw(base_size=20) + theme_grey(20)
plot(g)
buf2 = paste("n",m_n,"s",m_sig_GP*10,"_MLE.eps", sep="")
dev.copy2eps(file=buf2, width=8, height=6)

# PMLE
x <- data.frame(
  x = as.factor(c(rep("MLE",iter.max),rep("GCV",iter.max),
                  rep("mAIC",iter.max),rep("GIC",iter.max),rep("GBIC",iter.max))),
  y = c(resMSE[,1],resMSE[,7],resMSE[,8],resMSE[,9],resMSE[,10]) / (m_n*m_ni),
  z=1:500
)
g <- ggplot(x,aes(x=reorder(x,z),y=y), xlab="")
g <- g + geom_boxplot() + xlab("") + ylab("ASE") +
  theme_bw(base_size=20) + theme_grey(20)
plot(g)
buf2 = paste("n",m_n,"s",m_sig_GP*10,"_IC.eps", sep="")

