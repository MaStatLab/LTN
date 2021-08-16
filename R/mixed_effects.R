#' Gibbs sampler for cross-group comparison
# N -- number of samples
# p -- number of internal nodes
# g -- number of levels of the random effect. e.g., if the random effect is individual, then g is number of individuals
# r -- pre-specified number of factors. unwarranted feature.
#' @param Y -- (numeric) N*p, pre-order
# Xtest -- N*q1 design [matrix] of the covariate to test, q1 should be [the number of groups-1]. If do not -1, will have 2 beta's for a two-group comparison, and need to test the difference!
# Xadjust -- N*q2 design [matrix] of the covariate to adjust, should include a column of 1's
# refflabel -- vector of random effect labels with length g. For example, patient ID.
# c0,d0,c1,d1,nu -- hyperparameters
# niter -- number of Gibbs iterations
# SEED -- random seed for initialization
# save_alpha_only -- T/F, whether to only save alpha. If T, will save huge amount of memory.
# gprior_m -- Var(beta2)=mN(X^TX)^-1, gprior_m=1 implies unit information prior
# OUTPUT:
# **All in the same order as nodes**
# ALPHA -- list, samples of r*p factor loading matrices
# BETA1 -- list, samples of q1*p regression coefficients for testing
# BETA2 -- list, samples of q2*p regression coefficients for adjusting
# Z -- list, samples of [N*r] latent factors (in each iteration it's r*N)
# GAM -- list, samples of [p*g] random effects
# PSI -- list, samples of [N*p] psi's (in each iteration it's p*N)
# PHI_EPS -- p*niter, samples of phi_epsilon's
# PHI_GAM -- p*niter, samples of phi_gamma's
# PHI_PI -- q1*niter, samples of phi_pi's
# ACC1, ACC2 -- vectors of whether the new samples are accepted in MH steps
# PI1 -- list, samples of q1*p inclusion probabilities of beta1
# OMEGA2
# LAM2
# HYPERPARAMETERS:
# aj, mj
# t0=1, u0=0.001 -- $\phi_{beta_j}=\sim Ga(t,u)$, denoted as phi_pi in code
# reffcov -- 1: diagonal, 2: sparse inverse covariance
# a1,a2 -- hyperparameter in prior on factor loadings.
#           Dunson (2011) requires a2>b2+1 where delta2~Ga(a2,b2), and a2>3,a1>2
#           if a1a2=='hp', these values of a1,a2 here will not be used.
# a1a2 -- hp: use Ga(2,1) hyperprior on a1,a2; fix: fix a1, a2
#' @export
gibbs_crossgroup = function(N,
                            p,
                            g = NULL,
                            r,
                            YL,
                            Y,
                            Xtest,
                            Xadjust = NULL,
                            grouplabel = NULL,
                            c0 = 1,
                            d0 = 0.001,
                            c1 = 1,
                            d1 = 0.001,
                            nu = 3,
                            niter,
                            adjust = F,
                            reff = F,
                            reffcov = 1,
                            SEED = 1,
                            save_alpha_only = F,
                            gprior_m = 1,
                            pnull = 0.5,
                            a1=3,
                            a2=4,
                            a1a2='hp',
                            lambda_fixed) {
  # initialization
  set.seed(SEED)
  Y = t(Y)
  YL = t(YL)
  if (!adjust) {
    Xadjust = matrix(rep(1, N), ncol = 1)
  }
  X = cbind(Xtest, Xadjust)
  phi_eps = rep(1, p)
  PHI_EPS = phi_eps
  if (reff) {
    if (reffcov == 1) {
      phi_gam = rep(1, p)
      PHI_GAM = phi_gam
      OMEGA2 = NULL
    } else{
      PHI_GAM = NULL
      omega2 = stats::rWishart(1, p + 2, diag(p))[, , 1]
      OMEGA2 = list(omega2)
      if (lambda_fixed>0){
        lam2 = lambda_fixed
      }else{
        lam2=1
      }
      LAM2 = lam2
      tau2 = matrix(0, nrow = p, ncol = p)
      for (l in 2:p) {
        for (k in 1:(l - 1)) {
          mu_prime = abs(lam2/ omega2[k, l])
          ukl = statmod::rinvgauss(1, mu_prime, lam2 ^ 2)
          tau2[k, l] = 1 / ukl
          tau2[l, k] = 1 / ukl
        }
      }
    }
  } else {
    OMEGA2=NULL
    LAM2=NULL
    PHI_GAM = NULL
    grouplabel = rep(1, N)
  }
  psi = psi_mle(YL, Y) # p*N
  kappa = YL - Y / 2
  PSI = list(psi)
  wpg = apply(Y, 1:2, function(x) {
    BayesLogit::rpg(1, as.numeric(x + 0.001), 0)
  })
  # p*N, do not really need to initialize
  # rpg function cannot be applied to "integer" type parameters!!!
  # h>0
  WPG = list(wpg)
  if (r > 0) {
    z = matrix(stats::rnorm(N * r), r, N)
    Z = list(z)
  } else{
    z = matrix(0, 2, N)
    Z = list(z)
  }
  if (reff) {
    gam = matrix(stats::rnorm(p * g), p, g)
  } else {
    g = 1
    gam = matrix(0, p, g)
  }
  GAM = list(gam)
  A1=NULL
  A2=NULL
  if (r > 0) {
    alpha = matrix(stats::rnorm(r * p), r, p)
    ALPHA = list(alpha)
    delta = stats::rgamma(r, 2, 1)
    tau = cumprod(delta)
    TAU = tau
    phi = matrix(stats::rgamma(r * p, nu / 2, nu / 2), r, p) # r*p
    PHI = list(phi)
    if (a1a2=='hp'){
      a1 = stats::rgamma(1,2, 1)
      a2 = stats::rgamma(1,2, 1)
    }
  } else{
    A1=NULL
    A2=NULL
    TAU=NULL
    PHI=NULL
    alpha = matrix(0, 2, p)
    ALPHA = list(alpha)
  }
  if (reff) {
    nl = sapply(1:g, function(x) {
      sum(grouplabel == x)
    })
  }
  q = ncol(X)
  q1 = ncol(Xtest)
  q2 = ncol(Xadjust)
  beta = solve(t(X) %*% X,t(X) %*% t(psi))
  beta1 = matrix(as.vector(beta[1:q1, ]), nrow = q1)
  BETA1 = list(beta1)
  beta2 = matrix(as.vector(beta[(q1 + 1):(q1 + q2), ]), nrow = q2)
  BETA2 = list(beta2)
  aj = 1
  mj = 1 - (pnull) ^ (1 / p)
  pi1 = matrix(stats::rbeta(q1 * p,aj*mj,aj*(1-mj)), q1, p)
  PI1 = list(pi1) # q1*p inclusion probability of beta1
  t0 = 1
  u0 = 0.001
  phi_pi = stats::rgamma(q1, t0, u0)
  PHI_PI = phi_pi
  # sampling
  ACC1 = rep(0, niter)
  ACC2 = rep(0, niter)
  for (it in 2:niter) {
    #print(it)
    if (reff) {
      # update gam
      if (reffcov == 1) {
        omega2 = diag(phi_gam)
      }
      Clist = lapply(nl, function(x) {
        chol2inv(chol(x * diag(phi_eps, nrow = p) + omega2))
      }) # list of covariance matrices
      diff = psi - t(alpha) %*% z - t(X %*% beta) # the matrix is p*N, each column is a sample
      groupdiff = sapply(1:g, function(x) {
        rowSums(as.matrix(diff[, grouplabel == x]))
      }) # p*g
      gamlist = lapply(1:g, function(l) {
        m = Clist[[l]] %*% diag(phi_eps, nrow = p) %*% matrix(groupdiff[, l], ncol =
                                                                1)
        return(t(mvtnorm::rmvnorm(1, m, Clist[[l]])))
      })
      gam = do.call(cbind, gamlist) # p*g
      if (!save_alpha_only) {
        GAM[[it]] = gam
      }
    }
    # update phi_eps (parallel)
    err = psi - gam[, grouplabel] - t(alpha) %*% z - t(X %*% beta) # p*N, each col is a sample
    ess=apply(err,1,function(x){
      x=x[x!=0]
      lx2=2*log(abs(x))
      return(exp(logsumexp(lx2)-log(2)))
    }) # p-vec
    phi_eps = sapply(ess, function(x) {
      stats::rgamma(1, c0 + N / 2, d0 + x)
    })
    phi_eps[phi_eps < 10 ^ (-10)] = 10 ^ (-10)
    PHI_EPS = cbind(PHI_EPS, phi_eps)
    if (reff) {
      # update cov(gamma)
      if (reffcov == 1) {
        #diagonal
        # update phi_gam
        ss=apply(gam,1,function(x){
          x=x[x!=0]
          lx2=2*log(abs(x))
          return(exp(logsumexp(lx2)-log(2)))
        }) # p-vec
        phi_gam = sapply(ss, function(x) {
          stats::rgamma(1, c1 + g / 2, d1 + x)
        })
        PHI_GAM = cbind(PHI_GAM, phi_gam)
      } else{
        # sparse
        S2 = gam %*% t(gam)
        omegatau = update_omega_tau(S2, lam2, g, omega2, tau2)
        omega2 = omegatau[[1]]
        tau2 = omegatau[[2]]
        if (lambda_fixed<=0){
          lam2 = stats::rgamma(1, 1 + p * (p + 1) / 2, 0.01 + exp(logsumexp(log(abs(omega2)))-log(2)))
        }
        OMEGA2[[it]] = omega2
        LAM2[[it]] = lam2
      }
    }
    # update psi
    C = 1 / apply(wpg, 2, function(x) {
      x + phi_eps
    }) # p*N
    b = kappa + apply(t(alpha) %*% z + gam[, grouplabel] + t(X %*% beta), 2, function(x) {
      x * phi_eps
    }) # p*N
    m = C * b
    mc = rbind(m, C) # p*N
    psi = apply(mc, 2, function(x) {
      mvtnorm::rmvnorm(1, x[1:p], diag(x[-(1:p)]))
    }) # p*N
    if (!save_alpha_only) {
      PSI[[it]] = t(psi)
    }
    # update wpg
    wpg = mapply(rpg0, Y, psi) # vec(wpg)
    wpg = matrix(wpg, ncol = N, nrow = p)
    if (!save_alpha_only) {
      WPG[[it]] = wpg
    }
    if (r > 0) {
      # update z
      C = chol2inv(chol(alpha %*% diag(phi_eps) %*% t(alpha) + diag(r))) # universal C, unavoidable matrix inversion. Use chol to ensure symmetry.
      m = C %*% alpha %*% (diag(phi_eps) %*% (psi - gam[, grouplabel] -
                                                t(X %*% beta))) # r*N
      z = apply(m, 2, function(x) {
        mvtnorm::rmvnorm(1, x, C)
      }) #r*N
      z = matrix(z, r, N)
      if (!save_alpha_only) {
        Z[[it]] = t(z)
      }
      # update alpha
      alpha = matrix(0, r, p)
      for (j in 1:p) {
        D = diag(phi[, j] * tau, nrow = r)
        C = chol2inv(chol(D + phi_eps[j] * z %*% t(z))) # unavoidable matrix inversion
        m = phi_eps[j] * C %*% z %*% matrix(psi[j, ] - gam[j, grouplabel] -
                                              X %*% beta[, j], ncol = 1)
        alpha[, j] = mvtnorm::rmvnorm(1, m, C)
      }
      if (!save_alpha_only) {
        ALPHA[[it]] = alpha
      }
      # update phi (in the prior of alpha)
      gampara2 = nu/2 + exp(sweep(log(abs(alpha))*2,1,log(tau),"+"))/2 # r*p
      phi = apply(gampara2, 1:2, function(x) {
        stats::rgamma(1, (nu + 1) / 2, x)
      })
      if (!save_alpha_only) {
        PHI[[it]] = phi
      }
      # update tau
      # update delta1
      # gampara2 = sum(c(1, cumprod(delta[-1])) * rowSums(phi * alpha ^ 2)) /2 + 1
      gampara2=exp(logsumexp(c(0, cumsum(log(delta[-1])))+ apply(log(phi)+2*log(abs(alpha)),1,logsumexp)-log(2)))+1
      delta[1] = stats::rgamma(1, a1 + p * r / 2, gampara2)
      # update delta 2-r
      if (r > 1) {
        for (h in 2:r) {
          taulh = prod(delta[1:(h - 1)])
          # s = taulh * sum(phi[h, ] * alpha[h, ] ^ 2)
          s = taulh * exp(logsumexp(log(phi[h, ]) +2*log(abs(alpha[h, ]))))
          if (r > h) {
            for (l in (h + 1):r) {
              taulh = taulh * delta[l]
              # s = s + taulh * sum(phi[l, ] * alpha[l, ] ^ 2)
              s = s + taulh * exp(logsumexp(log(phi[l, ]) +2*log(abs(alpha[l, ]))))
            }
          }
          gampara2 = 1 + s / 2
          delta[h] = stats::rgamma(1, a2 + p * (r - h + 1) / 2, gampara2)
        }
      }
      # tau = cumprod(delta)
      tau = exp(cumsum(log(delta)))
      if (!save_alpha_only) {
        TAU = cbind(TAU, tau)
      }
      # update a1 with MH, proposal distr: prior
      if (a1a2=='hp'){
        at = stats::rgamma(1, 2, 1) # choice of hyperparameter is from Dunson 2011
        # acc = min(1, gamma(a1) / gamma(at) * delta[1] ^ (at - a1))
        acc=exp(min(0,lgamma(a1)-lgamma(at)+(at-a1)*log(delta[1])))
        u = stats::runif(1)
        if (u < acc) {
          a1 = at
          if (!save_alpha_only) {
            ACC1[it] = 1
          }
        }
        A1=c(A1,a1)
        # update a2 with MH
        at = stats::rgamma(1, 2, 1)
        acc = exp(min(0, (sum(log(delta[-1])) * at - (r - 1) * lgamma(at)) - (sum(log(delta[-1])) * a2 - (r - 1) * lgamma(a2))))
        u = stats::runif(1)
        if (u < acc) {
          a2 = at
          if (!save_alpha_only) {
            ACC2[it] = 1
          }
        }
        A2=c(A2,a2)
      }
    }
    # update beta1 then pi1
    for (k in 1:p) {
      for (j in 1:q1) {
        # s1 = sqrt(1 / (phi_pi[j] + sum(Xtest[, j] ^ 2) * phi_eps[k]))
        N1=sum(Xtest[, j] )
        ls1=-logsumexp(c(log(phi_pi[j]),log(N1)+log(phi_eps[k])))/2
        s1=exp(ls1)
        if (q1 > 1) {
          b1 = phi_eps[k] * sum(
            Xtest[, j] * (
              t(psi)[, k] - t(gam[, grouplabel])[, k] - as.matrix(as.vector(Xtest[, -j]), ncol =
                                                                    q1 - 1) %*% as.matrix(as.vector(beta1[-j, k]), nrow = q1 - 1) - t(z) %*% alpha[, k] -
                Xadjust %*% matrix(as.vector(beta2[, k]), ncol = 1)
            )
          )
        } else {
          b1 = phi_eps[k] * sum(Xtest[, j] * (
            t(psi)[, k] - t(gam[, grouplabel])[, k] - t(z) %*% alpha[, k] - Xadjust %*%
              matrix(as.vector(beta2[, k]), ncol = 1)
          ))
        }
        # prob0 = (1 - pi1[j, k]) / (1 - pi1[j, k] + pi1[j, k] * s1 * sqrt(phi_pi[j]) * exp(b1 ^ 2 * s1 ^ 2 / 2))
        lprob0=log(1-pi1[j, k])-logsumexp(c(log(1-pi1[j, k]),log(pi1[j, k])+ls1+log(phi_pi[j])/2+b1^2*s1^2/2))
        prob0=exp(lprob0)
        if (is.infinite(logsumexp(c(log(1-pi1[j, k]),log(pi1[j, k])+ls1+log(phi_pi[j])/2+b1^2*s1^2/2)))){
          warning('infinite logsumexp in prob0')
        }
        i0 = stats::rbinom(1, 1, prob0) # indicator of 0
        if (i0 == 1) {
          beta1[j, k] = 0
          pi1[j, k] = 0
          pi1[j, k] = stats::rbeta(1, aj * mj, aj * (1 - mj) + 1)
        } else {
          beta1[j, k] = stats::rnorm(1, b1 * s1 ^ 2, s1) # stats::rnorm(n,mean,SD), mvtnorm::rmvnorm(n,mean,var)
          pi1[j, k] = stats::rbeta(1, aj * mj + 1, aj * (1 - mj))
        }
      }
    }
    BETA1[[it]] = beta1
    PI1[[it]] = pi1
    # update beta2
    invXTX = chol2inv(chol(t(Xadjust) %*% Xadjust)) # unavoidable matrix inversion
    m_common = invXTX %*% t(Xadjust) %*% (t(psi) -t(z) %*% alpha - t(gam[, grouplabel]) - Xtest %*% beta1)
    m = m_common %*% diag(phi_eps / (phi_eps + 1 / N / gprior_m)) # q2*p
    beta2 = t(mvtnorm::rmvnorm(p, rep(0, q2), invXTX)) %*% diag(1 / sqrt(phi_eps +
                                                                  1 / N / gprior_m)) + m # q2*p
    beta2 = matrix(as.vector(beta2), nrow = q2)
    if (!save_alpha_only) {
      BETA2[[it]] = beta2
    }
    beta = rbind(beta1, beta2)
    # update phi_pi_j's
    nj_prime = rowSums(beta1 != 0)
    # phi_pi = mapply(function(x, y) { stats::rgamma(1, x, y) }, t0 + nj_prime / 2, u0 + rowSums(beta1 ^ 2) / 2)
    sumbeta2=apply(beta1,1,function(x){
      if (all(x==0)){
        return(0)
      } else{
        return(exp(logsumexp(2*log(abs(x[x!=0])))-log(2)))
      }
    })
    phi_pi = mapply(function(x, y) { stats::rgamma(1, x, y) }, t0 + nj_prime / 2, u0 +  sumbeta2)
    PHI_PI = cbind(PHI_PI, phi_pi)
  }
  if (!save_alpha_only) {
    return(
      list(
        ALPHA = ALPHA,
        BETA1 = BETA1,
        BETA2 = BETA2,
        Z = Z,
        GAM = GAM,
        PSI = PSI,
        PHI_EPS = PHI_EPS,
        PHI_GAM = PHI_GAM,
        OMEGA2 = OMEGA2,
        ACC1 = ACC1,
        ACC2 = ACC2,
        PHI_PI = PHI_PI,
        PI1 = PI1,
        LAM2 = LAM2,
        A1=A1,
        A2=A2,
        TAU=TAU,
        PHI=PHI
      )
    )
  } else{
    return(list(BETA1 = BETA1))
  }
}

# Calculate MLE of psi for initializing the gibbs sampler
# INPUT:
# yl,y -- vectors, single observation
# OUTPUT:
# psi -- (vector) MLE of psi in the same order as y
psi_mle = function(yl, y) {
  ratio = yl / y
  ratio[is.na(ratio)] = 0.5
  ratio[ratio == 0] = 10 ^ (-5)
  ratio[ratio == 1] = 1 - 10 ^ (-5)
  psi = stats::qlogis(ratio)
  return(psi)
}
# a function called in gibbs sampler
rpg0 = function(h, z) {
  h = as.numeric(h)
  z = as.numeric(z)
  if (h == 0) {
    x = 0
  }
  else{
    x = BayesLogit::rpg(1, h, z)
  }
  while (is.na(x)) {
    if (h == 0) {
      x = 0
    }
    else{
      x = BayesLogit::rpg.gamma(1, h, z)
    }
  }
  return(x)
}

# update omega with BGLA prior
update_omega_tau = function(S, lambda, nsim, omega, tau) {
  p = nrow(S)
  for (k in 1:p) {
    #(b)
    s22 = S[k, k]
    gam = stats::rgamma(1, nsim / 2 + 1, (s22 + lambda) / 2)
    Dtau = diag(tau[-k, k])
    omega11 = omega[-k, -k]
    C = chol2inv(chol((s22 + lambda) * chol2inv(chol(omega11)) + chol2inv(chol(Dtau))))
    s21 = S[-k, k]
    beta_mean = -C %*% s21
    be = mvtnorm::rmvnorm(1, beta_mean, C)
    #(c)
    omega[-k, k] = be
    omega[k, -k] = t(be)
    omega[k, k] = gam + be %*% chol2inv(chol(omega[-k, -k])) %*% t(be)
  }
  # step2
  for (l in 2:p) {
    for (k in 1:(l - 1)) {
      mu_prime = sqrt(lambda ^ 2 / omega[k, l] ^ 2)
      ukl = statmod::rinvgauss(1, mu_prime, lambda ^ 2)
      tau[k, l] = 1 / ukl
      tau[l, k] = 1 / ukl
    }
  }
  return(list(omega, tau))
}

# log(sum(exp(x)))
logsumexp=function(x){ # x is a vector
  M=max(x)
  return(M+log(sum(exp(x-M))))
}
