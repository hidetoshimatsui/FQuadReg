###' GPR estimation (MLE)
#' [in] Z: predictor (matrix)
#' [in] y: response (matrix)
#' [in] Psi: Basis functions (matrix)
CalcGPR = function(Z, y, Psi)
{
  #initial values
  nTheta = matrix(rnorm(m_my*ncol(Z)), nr=m_my, nc=ncol(Z))
  Theta = nTheta+1
  nnu = c(1,1,1)*1e-2
  nu = nnu+1
  Kmat = matrix(nr=m_ni, nc=m_ni)
  dSmat = array(0, dim=c(m_ni, m_ni, 3))  #1st derivative of Sigma
  ddSmat = array(0, dim=c(m_ni, m_ni, 3, 3))  #2nd derivative of Sigma
  NN=0
  
  while(max(abs(Theta-nTheta)) > 1e-4){
    NN=NN+1
    Theta = nTheta

    #Create matrix K
    for(j in 1:m_ni){
      for(k in j:m_ni){
        Kmat[j, k] = Kmat[k, j] =  nnu[1]*exp(-0.5*nnu[2]*(m_ti[j]-m_ti[k])^2)
      }
    }
    Sigma = Kmat + diag(nnu[3], m_ni) #Sigma

    # Update Theta
    buf0 = buf1 = 0
    for(i in 1:m_n){
      ni = m_t[m_t1[i]:(m_t1[i+1]-1)] #time points for i-th subject
      buf01 = t(Z[i, ]) %x% Psi[ni, ] #Kronecker prod
      buf02 = solve(Sigma[ni, ni])        #Inverse of covariance matrix
      buf0 = buf0 + t(buf01) %*% buf02 %*% buf01  #(XX)
      buf1 = buf1 + t(buf01) %*% buf02 %*% y[i, ] #Xy
      #print(c(i, det(buf0)))
    }
    nvecTheta= ginv(buf0) %*% buf1
    nTheta = matrix(nvecTheta, nr=m_my)

    # Update nu
    nu = nnu+1
    while(max(abs(nu-nnu))>1e-4){
      nu = nnu
      
      #0th and 1st derivative
      for(j in 1:m_ni){
        for(k in j:m_ni){
          buf1 = (m_ti[j]-m_ti[k])^2
          buf2 = exp(-0.5*nu[2]*buf1)
          Kmat[j, k] = Kmat[k, j] =  nu[1]*buf2
          dSmat[j, k, 1] = dSmat[k, j, 1] = buf2
          dSmat[j, k, 2] = dSmat[k, j, 2] = -0.5*nu[1] * buf1 * buf2
          ddSmat[j, k, 1, 2] = ddSmat[k, j, 2, 1] = -0.5 * buf1 * buf2  
          ddSmat[j, k, 2, 2] = 0.25 * nu[1] * buf1^2 * buf2
          #Other elements are zero
        }
      }
      Sigma = Kmat + diag(nu[3], m_ni) #Sigma
      dSmat[, , 3] = diag(1, m_ni)     #Sigma(1st derivative)

      dl = rep(0, 3)
      ddl = matrix(0, nr=3, nc=3)
      for(i in 1:m_n){
        ni = m_t[m_t1[i]:(m_t1[i+1]-1)] #time points for i-th subject
        Siginvi = solve(Sigma[ni, ni])
        
        r = y[i, ] - Psi[ni, ] %*% nTheta %*% Z[i, ]
        alpi = Siginvi %*% r
        alpi2buf = alpi%*%t(alpi)
        for(l in 1:3){
          #1st derivative of likelihood
          dl[l] = dl[l] + 0.5 * (sum(diag((alpi2buf - Siginvi) %*% dSmat[ni,ni,l])))
          
          #2nd derivative of likelihood
          for(ll in l:3){
            Amat = dSmat[ni, ni, l] %*% Siginvi %*% dSmat[ni, ni, ll]
            ddl[l, ll] = ddl[l, ll] + 
              0.5 * sum(diag((alpi2buf - Siginvi) %*% (ddSmat[,,l,ll] - Amat)- 
                              alpi2buf %*% Amat))
          }
        }
      }
      for(l in 1:3){
        for(ll in l:3){
          ddl[ll, l] = ddl[l, ll]
        }
      }
      #Newton method
      nnu = nu - solve(ddl)%*%dl
      #print(nnu[3])
    } # end while nu
  #print(NN)
  }
  return(list(Theta=nTheta, nu = nnu))
}

###' GPR estimation (Penalized MLE)
#' [in] Z: predictor (matrix)
#' [in] y: response (matrix)
#' [in] Psi: Basis functions (matrix)
#' [in] lambda: regularization parameter
CalcPGPR = function(Z, y, Psi, lambda)
{
  #initial values
  nTheta = matrix(rnorm(m_my*ncol(Z)), nr=m_my, nc=ncol(Z))
  Theta = nTheta+1
  nnu = c(1,1,1)*1e-2
  nu = nnu+1
  Kmat = matrix(nr=m_ni, nc=m_ni)
  dSmat = array(0, dim=c(m_ni, m_ni, 3))  #1st derivative of Sigma
  ddSmat = array(0, dim=c(m_ni, m_ni, 3, 3))  #2nd derivative of Sigma
  NN=0

  #regularization term
  Omegax = Penmat(m_mx) 
  Omegay = Penmat(m_my)
  mcol = 1 + m_mx + m_mx^2
  pen1 = matrix(0, nr = mcol, nc = mcol)
  pen1[2:(m_mx+1), 2:(m_mx+1)] = Omegax
  pen1[(m_mx+2):nrow(pen1), (m_mx+2):nrow(pen1)] = 
    Omegax%x%diag(1, m_mx) + diag(1, m_mx)%x%Omegax
  pen = diag(1, m_my)%x%pen1 + Omegay%x%diag(1, mcol)

  while(max(abs(Theta-nTheta)) > 1e-4 && NN<100){
    NN=NN+1
    Theta = nTheta
    
    #Create matrix K
    for(j in 1:m_ni){
      for(k in j:m_ni){
        Kmat[j, k] = Kmat[k, j] =  nnu[1]*exp(-0.5*nnu[2]*(m_ti[j]-m_ti[k])^2)
      }
    }
    Sigma = Kmat + diag(nnu[3], m_ni) #Sigma
    
    # Update Theta
    buf0 = buf1 = 0
    for(i in 1:m_n){
      ni = m_t[m_t1[i]:(m_t1[i+1]-1)] #time points for i-th subject
      zPsi = t(Z[i, ]) %x% Psi[ni, ]  #Kronecker prod
      iSig = solve(Sigma[ni, ni])        #inverse of covariance matrix
      buf0 = buf0 + t(zPsi) %*% iSig %*% zPsi  #(XX)
      buf1 = buf1 + t(zPsi) %*% iSig %*% as.numeric(y[i, ]) #Xy
      #print(c(i, det(buf0)))
    }
    nvecTheta= solve(buf0 + lambda*nrow(Z)*pen) %*% buf1
    nTheta = matrix(nvecTheta, nr=m_my)
    
    # Update nu
    #browser()
    nu = nnu+1
    while(max(abs(nu-nnu))>1e-4){
      nu = nnu
      
      #0th and 1st derivative
      for(j in 1:m_ni){
        for(k in j:m_ni){
          buf1 = (m_ti[j]-m_ti[k])^2
          buf2 = exp(-0.5*nu[2]*buf1)
          Kmat[j, k] = Kmat[k, j] =  nu[1]*buf2
          dSmat[j, k, 1] = dSmat[k, j, 1] = buf2
          dSmat[j, k, 2] = dSmat[k, j, 2] = -0.5*nu[1] * buf1 * buf2
          ddSmat[j, k, 1, 2] = ddSmat[k, j, 2, 1] = -0.5 * buf1 * buf2  
          ddSmat[j, k, 2, 2] = 0.25 * nu[1] * buf1^2 * buf2
          #Other elements are zero
        }
      }
      Sigma = Kmat + diag(nu[3], m_ni) #Sigma
      dSmat[, , 3] = diag(1, m_ni)     #Sigma (1st derivative)
      
      dl = rep(0, 3)
      ddl = matrix(0, nr=3, nc=3)
      for(i in 1:m_n){
        ni = m_t[m_t1[i]:(m_t1[i+1]-1)] #time points for i-th subject
        Siginvi = solve(Sigma[ni, ni])
        
        r = as.numeric(y[i, ]) - Psi[ni, ] %*% nTheta %*% Z[i, ]
        alpi = Siginvi %*% r
        alpi2buf = alpi%*%t(alpi)
        for(l in 1:3){
          #1st derivative of likelihood
          dl[l] = dl[l] + 0.5 * (sum(diag((alpi2buf - Siginvi) %*% dSmat[ni,ni,l])))
          
          #2nd derivative of likelihood
          for(ll in l:3){
            Amat = dSmat[ni, ni, l] %*% Siginvi %*% dSmat[ni, ni, ll]
            ddl[l, ll] = ddl[l, ll] + 
              0.5 * sum(diag((alpi2buf - Siginvi) %*% (ddSmat[,,l,ll] - Amat)- 
                               alpi2buf %*% Amat))
          }
        }
      }
      for(l in 1:3){
        for(ll in l:3){
          ddl[ll, l] = ddl[l, ll]
        }
      }
      #Newton method
      nnu = nu - solve(ddl)%*%dl
      #print(nnu[3])
    } # end while nu
    #print(NN)
  }
  return(list(Theta=nTheta, nu = nnu))
}


###' AIC and BIC for GPR
#' [in] Z: predictor (matrix)
#' [in] y: response (matrix)
#' [in] Psi: Basis functions (matrix)
#' [in] est: Estimator
IC_GPR = function(Z, y, Psi, est){
  
  loglik = 0
  Sigma = matrix(nr=m_ni, nc=m_ni)
  for(i in 1:m_ni){
    for(j in 1:m_ni){
      Sigma[i,j] = est$nu[1] * exp(-0.5*est$nu[2]*(m_ti[i]-m_ti[j])^2)
    }
  }
  Sigma = Sigma + diag(est$nu[3], m_ni) 
  
  #Log-likelihood
  for(i in 1:m_n){
    ni = m_t[m_t1[i]:(m_t1[i+1]-1)] #time points for i-th subject
    r = y[i, ] - Psi[ni, ] %*% est$Theta %*% Z[i, ]
    r = as.numeric(r)
    loglik = loglik - 0.5*m_ni*log(2*pi) - 0.5*log(det(Sigma[ni, ni])) -
             0.5 * t(r) %*% solve(Sigma[ni, ni]) %*% r
  }
  loglik = as.numeric(loglik)
  
  #degrees of freedom
  df = prod(dim(est$Theta)) + length(est$nu)
  
  AIC = -2*loglik + 2*df
  BIC = -2*loglik + log(m_n)*df
  return(list(AIC=AIC, BIC=BIC))
}


###' GCV mAIC, mBIC, GIC for GPR
#' [in] Z: predictor (matrix)
#' [in] y: response (matrix)
#' [in] Psi: Basis functions (matrix)
#' [in] est: Estimator
IC_PGPR = function(Z, y, Psi, est, lambda){
  
  loglik = 0
  thetavec = as.vector(est$Theta)  #m_my x (1+m_mx+m_mx^2)
  nu = est$nu
  Sigma = matrix(nr=m_ni, nc=m_ni)
  Kmat = Sigma
  for(i in 1:m_ni){
    for(j in 1:m_ni){
      Sigma[i,j] = est$nu[1] * exp(-0.5*est$nu[2]*(m_ti[i]-m_ti[j])^2)
    }
  }
  Sigma = Sigma + diag(est$nu[3], m_ni) 
  
  #regularization term
  Omegax = Penmat(m_mx) 
  Omegay = Penmat(m_my)
  pen1 = matrix(0, nr = 1 + m_mx + m_mx^2, nc = 1 + m_mx + m_mx^2)
  pen1[2:(m_mx+1), 2:(m_mx+1)] = Omegax
  pen1[(m_mx+2):nrow(pen1), (m_mx+2):nrow(pen1)] = 
    Omegax%x%diag(1, m_mx) + diag(1, m_mx)%x%Omegax
  pen = diag(1, m_my)%x%pen1 + Omegay%x%diag(1, 1+m_mx+m_mx^2)

  #Log-likelihood
  for(i in 1:m_n){
    ni = m_t[m_t1[i]:(m_t1[i+1]-1)] #time points for i-th subject
    r = as.numeric(y[i, ]) - Psi[ni, ] %*% est$Theta %*% Z[i, ]
    r = as.numeric(r)
    loglik = loglik - 0.5*m_ni*log(2*pi) - 0.5*log(det(Sigma[ni, ni])) -
      0.5 * t(r) %*% solve(Sigma[ni, ni]) %*% r
  }
  loglik = as.numeric(loglik)
  
  #degrees of freedom
  buf0 = buf1 = 0
  for(i in 1:m_n){
    ni = m_t[m_t1[i]:(m_t1[i+1]-1)] #time points for i-th subject
    Xi = t(Z[i, ]) %x% Psi[ni, ] #Kronecker prod
    iSig = solve(Sigma[ni, ni])        #inverse of covariance matrix
    buf0 = buf0 + t(Xi) %*% iSig %*% Xi  #(XX)
    #print(c(i, det(buf0)))
  }
  #hat matrix
  Hatmat = solve(buf0 + lambda*nrow(Z)*pen) %*% buf0
  # degrees of freedom
  df = sum(diag(Hatmat)) + length(est$nu)

  # derivative of sigma w.r.t nu
  dSmat = array(0, dim=c(m_ni, m_ni, 3))  #1st derivative of Sigma
  ddSmat = array(0, dim=c(m_ni, m_ni, 3, 3))  #2nd derivative of Sigma
  #0th and 1st derivative
  for(j in 1:m_ni){
    for(k in j:m_ni){
      buf1 = (m_ti[j]-m_ti[k])^2
      buf2 = exp(-0.5*nu[2]*buf1)
      Kmat[j, k] = Kmat[k, j] =  nu[1]*buf2
      dSmat[j, k, 1] = dSmat[k, j, 1] = buf2
      dSmat[j, k, 2] = dSmat[k, j, 2] = -0.5*nu[1] * buf1 * buf2
      ddSmat[j, k, 1, 2] = ddSmat[k, j, 2, 1] = -0.5 * buf1 * buf2
      ddSmat[j, k, 2, 2] = 0.25 * nu[1] * buf1^2 * buf2
    }
  }
  Sigma = Kmat + diag(nu[3], m_ni) #Sigma
  dSmat[, , 3] = diag(1, m_ni)     #Sigma‚P‰ñ”÷•ª

  # 1st and 2nd derivative of penalized log-likelihood
  dl_Theta = rep(0, length(thetavec))
  R11 = array(0, dim=c(length(thetavec), length(thetavec)))
  R12 = array(0, dim=c(length(thetavec),3))
  Q11 = R11
  Q21 = t(R12)
  dl_nu = rep(0, 3)
  R22 = matrix(0, nr=3, nc=3)
  for(i in 1:m_n){
    ni = m_t[m_t1[i]:(m_t1[i+1]-1)] #time points for i-th subject
    Siginvi = solve(Sigma[ni, ni])
    r = as.numeric(y[i, ]) - Psi[ni, ] %*% est$Theta %*% Z[i, ]
    alpi = Siginvi %*% r  
    alpi2buf = alpi%*%t(alpi)  
    Xi = t(Z[i, ]) %x% Psi[ni, ]
    
    # 1st and 2nd derivative w.r.t. Theta
    dl_Theta = dl_Theta + t(Xi)%*%alpi - lambda*nrow(Z)*pen%*%thetavec
    R11 = R11 - t(Xi)%*%Siginvi%*%Xi - lambda*nrow(Z)*pen  #dll_Theta
    # 1st and 2nd derivative w.r.t. nu
    for(l in 1:3){
      #1st derivative of likelihood
      dl_nu[l] = dl_nu[l] + 0.5 * (sum(diag((alpi2buf - Siginvi) %*% dSmat[ni,ni,l])))

      #2nd derivative of likelihood
      for(ll in l:3){
        Amat = dSmat[ni, ni, l] %*% Siginvi %*% dSmat[ni, ni, ll]
        R22[l, ll] = R22[l, ll] +
          0.5 * sum(diag((alpi2buf - Siginvi) %*% (ddSmat[,,l,ll] - Amat)-
                           alpi2buf %*% Amat))
      }
    }
    #2nd derivative w.r.t. Theta and nu
    for(l in 1:3){
      R12[,l] = R12[,l] - t(Xi)%*%Siginvi%*%dSmat[,,l]%*%alpi
    }
    
    Q11 = Q11 + t(Xi)%*%alpi2buf%*%Xi - lambda*nrow(Z)*pen%*%thetavec%*%t(alpi)%*%Xi
    for(l in 1:3){
      Q21[l,] = Q21[l,] + 0.5 * (sum(diag((alpi2buf - Siginvi) %*% dSmat[ni,ni,l]))) * t(Xi)%*%alpi
    }
  }
  Q12 = dl_Theta%*%t(dl_nu)
  Q22 = dl_nu%*%t(dl_nu)
  for(l in 1:3){
    for(ll in l:3){
      R22[ll, l] = R22[l, ll]
    }
  }

  # R and Q matrix for GIC and GBIC (times n)
  Rmat = -rbind(cbind(R11, R12), cbind(t(R12), R22))
  Qmat =  rbind(cbind(Q11, Q12), cbind(Q21, Q22))
  degen = nrow(pen) - qr(pen)$rank

  #print(c(-2*loglik, 2*df))
  
  r = 0
  for(i in 1:m_n){
    ni = m_t[m_t1[i]:(m_t1[i+1]-1)] #time points for i-th subject
    r = r + (as.numeric(y[i, ]) - Psi[ni, ] %*% est$Theta %*% Z[i, ])^2
  }
  
  #Calculate IC
  GCV = sum(r) / (m_n*(1-df/m_n)^2)
  AIC = -2*loglik + 2*df

  GIC = -2*loglik + 2*sum(diag(solve(Rmat)%*%Qmat))
  GBIC = -2*loglik + nrow(Z)*lambda*t(thetavec)%*%pen%*%thetavec -
          (nrow(pen)-degen)*log(lambda) + degen*log(nrow(Z)) +
          sum(log(eigen(Rmat)$values))
  return(c(GCV, AIC, GIC, GBIC))
}







###' Generate GP
GP = function(n, r=1){
  mu = rep(0, n)
  gpmat = matrix(nr=n, nc=n)
  for(i in 1:n){
    for(j in 1:n){
      gpmat[i,j] = r * exp(-(i-j)^2/(2*10^2))
    }
  }
  gpval = mvrnorm(1, mu, gpmat)
  return(gpval)
}


###' Penalty matrix
Penmat = function(m){
  D = matrix(0, nrow=m-2, ncol=m)
  for(j in 1:(m-2)){
    D[j, j] = 1
    D[j, j+1] = -2
    D[j, j+2] = 1
  }
  return (t(D)%*%D)
}

