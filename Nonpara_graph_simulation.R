library(stats)
library(GUEST)
library(MASS)
library(glasso)
library(huge)
library(pROC)
library(clime)

i <- 1i  # imaginary unit

# --- Empirical characteristic function ---
emp_char <- function(x, u) {
  sapply(u, function(uu) mean(exp(i * uu * x)))
}

# --- Gaussian kernel Fourier factor ---
phi_kernel_normal <- function(u, h) exp(-0.5 * (h^2) * u^2)

# --- Deconvolution density estimator for one margin ---
deconv_adj_density <- function(Xstar, phi_eps, h = NULL, x_grid = NULL,
                               u_max = NULL, du = NULL, debug = FALSE) {
  n <- length(Xstar)
  if (is.null(h)) h <- bw.SJ(Xstar)
  if (h <= 0 || !is.finite(h)) h <- 1.06 * sd(Xstar) * n^(-1/5)
  
  if (is.null(u_max)) u_max <- 50 / max(h, 1e-6)
  if (is.null(du)) du <- u_max / 5000
  u <- seq(-u_max, u_max, by = du)
  
  # error characteristic function
  pe <- phi_eps(u)
  eps <- 1e-12
  too_small <- Mod(pe) < eps | is.na(pe)
  if (any(too_small)) pe[too_small] <- eps * exp(i * Arg(pe[too_small]))
  
  # empirical phi_{X*} with kernel
  phiXstar_hat <- phi_kernel_normal(u, h) * emp_char(Xstar, u)
  ratio <- phiXstar_hat / pe
  
  # x grid
  if (is.null(x_grid)) {
    rng <- range(Xstar)
    pad <- max(3*h, diff(rng)*0.2)
    x_grid <- seq(rng[1]-pad, rng[2]+pad, length.out = 1024)
  }
  
  # inverse Fourier
  inv_int <- function(x) {
    integrand <- Re(exp(-i*u*x) * ratio)
    integrand[is.na(integrand) | is.infinite(integrand)] <- 0
    du * (sum(integrand) - 0.5*(integrand[1]+integrand[length(integrand)])) / (2*pi)
  }
  f_vals <- vapply(x_grid, inv_int, numeric(1))
  f_vals[is.na(f_vals) | f_vals < 0] <- 0
  
  # normalize
  dx <- diff(x_grid)
  area <- sum((f_vals[-1] + f_vals[-length(f_vals)])*dx/2, na.rm = TRUE)
  if (!is.finite(area) || area <= 0) area <- 1
  f_vals <- f_vals / area
  
  # cumulative distribution function
  F_vals <- numeric(length(x_grid))
  F_vals[-1] <- cumsum((f_vals[-1]+f_vals[-length(f_vals)])*dx/2)
  F_vals[F_vals>1] <- 1
  F_interp <- approxfun(x_grid, F_vals, rule=2, ties="ordered")
  
  if (debug) {
    message("---- DEBUG ----")
    message("h=", signif(h,4), " u_max=", u_max, " dx=", signif(dx[1],4))
    message("Density range: [", signif(min(f_vals),4), ",", signif(max(f_vals),4), "]")
    message("Area before normalization=", signif(area,4))
    message("----------------")
  }
  
  list(x=x_grid, f=f_vals, F=F_vals, F_eval=F_interp, h=h, u=u)
}

# --- Safe distance correlation with regularization ---
distance_correlation <- function(X, Y, eps=1e-3) {
  n <- length(X)
  if (length(Y)!=n) stop("X and Y must be same length")
  
  A <- as.matrix(dist(X))
  B <- as.matrix(dist(Y))
  
  # double centering
  A <- A - rowMeans(A) - colMeans(A) + mean(A)
  B <- B - rowMeans(B) - colMeans(B) + mean(B)
  
  dcov2 <- mean(A*B)
  dvarX <- mean(A*A)
  dvarY <- mean(B*B)
  
  # regularization: prevent zero variance
dvarX <- max(dvarX, eps)
  dvarY <- max(dvarY, eps)
  
  ratio <- dcov2 / sqrt(dvarX*dvarY)
  
  # regularization: if ratio <= 0 (numerical issues), set to tiny positive
  if (!is.finite(ratio) || ratio <= 0) ratio <- eps
  
  sqrt(ratio)
}

# --- Multivariate X*: compute adjusted F and omega matrix ---
estimate_omega <- function(Xstar, phi_eps_list, h_vec=NULL, debug=FALSE) {
  Xstar <- as.matrix(Xstar)
  n <- nrow(Xstar)
  p <- ncol(Xstar)
  if (length(phi_eps_list)!=p) stop("phi_eps_list must have length p")
  
  Fadj <- matrix(NA_real_, n, p)
  deconv_list <- vector("list", p)
  
  for (j in 1:p) {
    h_j <- if (is.null(h_vec)) NULL else h_vec[j]
    obj <- deconv_adj_density(Xstar[,j], phi_eps_list[[j]], h=h_j, debug=debug)
    deconv_list[[j]] <- obj
    # add tiny jitter to avoid exact constants
    Fadj[,j] <- obj$F_eval(Xstar[,j]) + rnorm(n,0,1e-12)
  }
  
  omega <- diag(1, p)
  for (s in 1:(p-1)) {
    for (t in (s+1):p) {
      omega[s,t] <- omega[t,s] <- distance_correlation(Fadj[,s], Fadj[,t])
    }
  }
  
  list(Fadj=Fadj, omega=omega, deconv=deconv_list)
}

# --- Demo function ---
demo <- function(n=200, p=3, sigma_eps=c(0.3,0.5,0.2), debug=TRUE) {
  set.seed(123)
  Xtrue <- matrix(rnorm(n*p), n, p)
  Xstar <- Xtrue + matrix(rnorm(n*p, sd=sigma_eps), n, p)
  
  phi_list <- lapply(sigma_eps, function(s) function(u) exp(-0.5*s^2*u^2))
  
  fit <- estimate_omega(Xstar, phi_list, debug=debug)
  
  cat("Omega matrix:\n")
  print(round(fit$omega,3))
  
  # Plot first margin density
  plot(fit$deconv[[1]]$x, fit$deconv[[1]]$f, type="l", col="blue", lwd=2,
       main="Adjusted density of X1*", xlab="x", ylab="Density")
  lines(density(Xtrue[,1]), col="red", lty=2, lwd=2)
  lines(density(Xstar[,1]), col="green", lty=3, lwd=2)
  legend("topright", legend=c("Adjusted f_adj","True X1","Observed X1*"),
         col=c("blue","red","green"), lwd=2, lty=c(1,2,3))
  
  invisible(fit)
}


##################################
### Kendall's tau
##################################


Kendall = function(x) {

p = dim(Xstar)[2]
KM = matrix(0,p,p)
for(i in 1:p){
for(j in 1:p) {
KM[i,j] = sin( (pi/2) * cor.test(Xstar[,i],Xstar[,j],method="kendall")$estimate  )

}
}

return(KM)

}

Spearman = function(x) {

p = dim(Xstar)[2]
KM = matrix(0,p,p)
for(i in 1:p){
for(j in 1:p) {
KM[i,j] = 2*sin( (pi/6) * cor.test(Xstar[,i],Xstar[,j],method="spearman")$estimate  )

}
}

return(KM)

}


############################################################################
############################################################################
############################################################################
############################################################################


set.seed(123)
###
n <- 800
p <- 20
re = 50
PM = matrix(0,p,p); PM0 = matrix(0,p,p); PM_glasso = matrix(0,p,p)
PM_GUEST = matrix(0,p,p); PM_Kendall = matrix(0,p,p); PM_huge = matrix(0,p,p)
PM_Spearman = matrix(0,p,p); PM_clime = matrix(0,p,p)
for(ite in 1:re) {

# -----------------------
# Covariance generation for models 1-4
# -----------------------
generate_model1 <- function(p) {
  Theta <- matrix(0, p, p)
  diag(Theta) <- 1
  for(i in 1:(p-1)) {
    Theta[i, i+1] <- .5
    Theta[i+1, i] <- .5
  }
#  diag(Theta) = max(eigen(Theta)$value)
  solve(Theta)
}

generate_model2 <- function(p) {
  Theta <- matrix(0, p, p)
  diag(Theta) <- 1
  for(i in 1:(p-1)) Theta[i, i+1] <- Theta[i+1, i] <- 0.4
  for(i in 1:(p-2)) Theta[i, i+2] <- Theta[i+2, i] <- 0.2
  for(i in 1:(p-3)) Theta[i, i+3] <- Theta[i+3, i] <- 0.2
  solve(Theta)
}

generate_model3 <- function(p, num_hubs=16, links_per_hub=5) {
  Theta <- matrix(0, p, p)
  hub_nodes <- sample(1:p, num_hubs)
  for(h in hub_nodes) {
    neighbors <- sample(setdiff(1:p,h), links_per_hub)
    Theta[h, neighbors] <- 0.2
    Theta[neighbors, h] <- 0.2
  }
  diag(Theta) <- 1
  solve(Theta)
}

generate_model4 <- function(p) {
  Theta0 <- matrix(0, p, p)
  for(i in 1:(p-1)) {
    for(j in (i+1):p) {
      Theta0[i,j] <- Theta0[j,i] <- sample(c(0,0.2),1, prob=c(0.99,0.01))
    }
  }
  eig_min <- min(eigen(Theta0)$values)
  sigma_diag <- abs(eig_min) + 0.1
  Theta0 + diag(sigma_diag, p)
}




# -----------------------
# Gaussian simulation
# -----------------------
simulate_gaussian <- function(Sigma, n) {
  mvrnorm(n=n, mu=rep(0,ncol(Sigma)), Sigma=Sigma)
}

# -----------------------
# Safe Nonparanormal transformation
# -----------------------
simulate_nonparanormal_noNaN <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  X_trans <- matrix(0, n, p)
  
  for(j in 1:p) {
    xj <- X[,j]
    f_choice <- sample(1:5,1)
    
    if(f_choice == 1){          # identity
      X_trans[,j] <- xj
    } else if(f_choice == 2){   # log
      xj <- xj - min(xj) + 1e-3
      X_trans[,j] <- log(abs(xj))
    } else if(f_choice == 3){   # cube root
      X_trans[,j] <- sign(xj) * abs(xj)^(1/3)
    } else if(f_choice == 4){   # logit
      xj <- xj - min(xj)
      rng <- max(xj)
      if(rng == 0) rng <- 1
      xj <- xj / rng
      xj <- pmin(pmax(xj, 1e-6), 1-1e-6)
      X_trans[,j] <- log(abs(xj / (1-xj)))
    } else if(f_choice == 5){   # f5 safe
      xj <- xj - min(xj) + 1e-3
      X_trans[,j] <- ifelse(xj < 1, xj, log(abs(xj-1) + 1e-3) + 1)
    }
  }
  
  return(X_trans)
}

# -----------------------
# 生成 4 個模型的 Gaussian & Nonparanormal 資料
# -----------------------
Sigma1 <- generate_model1(p)


X1 <- simulate_gaussian(Sigma1, n)


X1b <- simulate_nonparanormal_noNaN(X1)


# -----------------------
# 檢查結果
# -----------------------
any(is.na(X1b))  # FALSE


dim(X1b)  # 300 x 100

#############################

root = function(r, L=100) {

X = seq(0,1,length=L)

f = function(x) {
return(
sqrt((x * asin(x) + sqrt(1-x^2) - x * asin(x/2) - sqrt(4-x^2) +1 ) /
(1+pi/3 - sqrt(3)) ) )
                           }

val = abs(f(X) - r)

return(  X[which(val == min(val))] )


}




sigma_eps = rep(1,p)
sigma_eps0 = rep(0,p)
Xstar <- X1b + matrix(rnorm(n*p, sd=sigma_eps), n, p)
  
  phi_list <- lapply(sigma_eps, function(s) function(u) exp(-0.5*s^2*u^2))
  phi_list0 <- lapply(sigma_eps0, function(s) function(u) exp(-0.5*s^2*u^2))  

  fit <- estimate_omega(Xstar, phi_list)
  fit0 <- estimate_omega(Xstar, phi_list0)

Rho = matrix(0,p,p); Rho0 = matrix(0,p,p)
Omega = fit$omega; Omega0 = fit0$omega
for(j in 1:p){
for(k in 1:p) {
Rho[j,k] = root(Omega[j,k])
Rho0[j,k] = root(Omega0[j,k])
}
}

PM = PM + abs(glasso(Rho,rho=0.1)$wi)
PM0 = PM0 + abs(glasso(Rho0,rho=0.15)$wi)
#PM = PM + abs(glasso(Omega,rho=0.1)$wi)
#PM0 = PM0 + abs(glasso(Omega0,rho=0.15)$wi)

PM_glasso = PM_glasso + glasso(cov(Xstar),rho=0.15)$wi
PM_GUEST = PM_GUEST + boost.graph(Xstar, ite1 = 20, thre = 0.1, sigma_e = sigma_eps[1],rep=1 )$w
PM_Kendall = PM_Kendall + glasso(Kendall(Xstar),rho=0.15)$wi
PM_Spearman = PM_Spearman + glasso(Spearman(Xstar),rho=0.15)$wi
#PM_huge = PM_huge + huge.glasso(cov(Xstar))$icov[[10]]
CLIME = clime(Xstar)
lambda0 = cv.clime(CLIME)$lambdaopt
lambda = which(CLIME$lambda == lambda0)
PM_clime = PM_clime + CLIME$Omega[[lambda]]






}
####################

  Theta <- matrix(0, p, p)
  diag(Theta) <- 1
  for(i in 1:(p-1)) {
    Theta[i, i+1] <- .2
    Theta[i+1, i] <- .2
  }

Theta[1:5,1:5]
pm=PM/re
pm[which(abs(pm)<0.1)]=0
pm0=PM0/re
pm0[which(abs(pm0)<0.15)]=0
pm_glasso = PM_glasso/re
pm_glasso[which(abs(pm_glasso)<0.15)]=0
pm_GUEST = PM_GUEST/re
pm_GUEST[which(abs(pm_GUEST)<0.15)]=0
pm_Kendall = PM_Kendall/re
pm_Kendall[which(abs(pm_Kendall)<0.15)]=0
pm_Spearman = PM_Spearman/re
pm_Spearman[which(abs(pm_Spearman)<0.15)]=0
#PM[which(PM<0.3)]=0
pm_clime = PM_clime/re
pm_clime[which(abs(pm_clime)<0.15)]=0

####################################3

max(abs(pm - Theta))
max(abs(pm0 - Theta))
max(abs(pm_glasso - Theta))
max(abs(pm_GUEST - Theta))
max(abs(pm_Kendall - Theta))
max(abs(pm_Spearman - Theta))


roc_PM = roc(Theta,PM/re)
roc_PM0 = roc(Theta,pm0)
roc_PM_glasso = roc(Theta,pm_glasso)
roc_PM_GUEST = roc(Theta,pm_GUEST)
roc_PM_Kendall = roc(Theta,pm_Kendall)
roc_PM_Spearman = roc(Theta,pm_Spearman)
roc_PM_clime = roc(Theta,pm_clime)
#roc_PM_Oracle = roc(Theta,Theta)

plot(roc_PM,col=1)
plot(roc_PM0,add = TRUE,col=2,lty=2)
plot(roc_PM_glasso,add = TRUE,col=3,lty=3)
plot(roc_PM_clime,add = TRUE,col=4,lty=4)
plot(roc_PM_Kendall,add = TRUE,col=5,lty=5)
plot(roc_PM_Spearman,add = TRUE,col=6,lty=6)
plot(roc_PM_GUEST,add = TRUE,col=7,lty=7)
#plot(roc_PM,col=1)

legend("bottomright",
legend = c(
"proposed","naive","glasso","GUEST",
"Kendall","Spearman","clime"),
col = c(1:7),lwd=2)

