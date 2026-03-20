####################
#### Data Setup ####
####################

data = read.table("C://GSE158588.csv",sep=",")
gene = data[,1]
id = data[1,]

data = data[,-1]; data = data[-1,]

data = matrix(as.numeric(unlist(data)),dim(data)[1],dim(data)[2])

sel = which(rowMeans(data)>=4000)

gene = gene[sel]
data = data[sel,]
dim(data)

################################
#### Methods Implementation ####
################################

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


Xstar <- t(data); Xstar = scale(Xstar)
n = dim(Xstar)[1] ; p = dim(Xstar)[2]   
sigma_eps = rep(.0001,p)
sigma_eps0 = rep(0,p)
  
  phi_list <- lapply(sigma_eps, function(s) function(u) exp(-0.5*s^2*u^2))
  phi_list0 <- lapply(sigma_eps0, function(s) function(u) exp(-0.5*s^2*u^2))  

  fit <- estimate_omega(Xstar, phi_list)
  fit0 <- estimate_omega(Xstar, phi_list0)

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



Rho = matrix(0,p,p); Rho0 = matrix(0,p,p)
Omega = fit$omega; Omega0 = fit0$omega
for(j in 1:p){
for(k in 1:p) {
Rho[j,k] = root(Omega[j,k])
Rho0[j,k] = root(Omega0[j,k])
}
}

PM = abs(glasso(Rho,rho=0.25)$wi)
PM0 =  abs(glasso(Rho0,rho=0.3)$wi)

PM_glasso =  glasso(cov(Xstar),rho=0.15)$wi
PM_GUEST =  boost.graph(Xstar, ite1 = 10, thre = 0.5, sigma_e = sigma_eps[1],rep=1 )$w
PM_GUEST_naive =  boost.graph(Xstar, ite1 = 10, thre = 0.5, sigma_e = sigma_eps0[1],rep=1 )$w
PM_Kendall =  glasso(Kendall(Xstar),rho=0.15)$wi
PM_Spearman =  glasso(Spearman(Xstar),rho=0.15)$wi
CLIME = clime(Xstar)
lambda0 = cv.clime(CLIME)$lambdaopt
lambda = which(CLIME$lambda == lambda0)
PM_clime = CLIME$Omega[[lambda]]


draw_graph = function(precision_matrix,label_name)
{
 net = precision_matrix
    net = network::network(net, directed = FALSE)
    network::network.vertex.names(net) = paste0("Y", network::network.vertex.names(net))
    graph = GGally::ggnet2(net, size = 10, node.color = "lightgray", 
        label = label_name, label.size = 3, mode = "circle")
    return(graph)
}

qq = 20
draw_graph(PM[1:qq,1:qq], gene[1:qq])
draw_graph(PM0[1:qq,1:qq], gene[1:qq])
draw_graph(PM_GUEST[1:qq,1:qq], gene[1:qq])
draw_graph(PM_glasso[1:qq,1:qq], gene[1:qq])
draw_graph(PM_Kendall[1:qq,1:qq], gene[1:qq])
draw_graph(PM_Spearman[1:qq,1:qq], gene[1:qq])
draw_graph(PM_clime[1:qq,1:qq], gene[1:qq])


(length(which(PM!=0)) - p)/2
(length(which(PM0!=0)) - p)/2
(length(which(PM_GUEST!=0)) - p)/2
(length(which(PM_glasso!=0)) - p)/2
(length(which(PM_Kendall!=0)) - p)/2
(length(which(PM_Spearman!=0)) - p)/2
(length(which(PM_clime!=0)) - p)/2

