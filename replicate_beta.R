rm(list = ls())

# the function is exactly same as gemma2::eigen2()
eigh <- function(spd, decreasing = FALSE) {
    foo <- eigen(spd, symmetric=TRUE)
    bar <- foo
    bar$values <- foo$values[order(foo$values, decreasing = decreasing)]
    bar$vectors <- foo$vectors[, order(foo$values, decreasing = decreasing)]
    return(bar)
}

##############################################################################
### load data ################################################################
##############################################################################

library(tidyverse)
library(data.table)

library(jjutil)
library(gemma2)

# Read Y
pheno <- readr::read_tsv(system.file("extdata", "mouse100.pheno.txt", package = "gemma2"), col_names = FALSE)
Y <- as.matrix(pheno[, c(1, 6)])
# Y <- as.matrix(pheno[, c(1, 2, 6)])

# Read Y
K <- readr::read_tsv(system.file("extdata", "mouse100.cXX.txt", package = "gemma2"), col_names = FALSE)[, 1:100]  %>% as.matrix

# Read X and single snp
geno <- readr::read_csv(system.file("extdata", "mouse100.geno.txt", package = "gemma2"), col_names = FALSE)
sel_idx <- which(geno$X1 %in% "rs8275764")
g1 <- geno[c(sel_idx), -c(1:3)] # first 3 columns are SNP annotations!
X <- t(as.matrix(g1))

# multiple snps
# sel_snps <- c("rs8275764", "rs6212654", "rs13477740")
# sel_idx <- match(sel_snps, geno$X1)
# g1 <- geno[c(sel_idx), -c(1:3)]
# X <- t(as.matrix(g1))

n_pheno <- ncol(Y)
n_indi <- nrow(Y)
n_snp <- ncol(X)

###########################################################################
### 1. data transformation based on the grm (i.e., kinship) ###############
###########################################################################
# we can find the reference for this transformation in the GEMMA paper

# eigendecomposition of K  
K_eigen <- eigh(K)
Kva <- K_eigen$values
Kve <- K_eigen$vectors
leftTransform <- Kve %>% t

# stadadization of Y
Y <- scale(Y)
Ystar <- leftTransform %*% Y


# X centering and transform
X <- scale(X, center = T, scale = F)
# We dont need to add the interaction term because we perform centeing the X matrix 
# X0_T <- leftTransform %*% matrix(rep(1, n_indi))
Xstar <- leftTransform %*% X

###########################################################################
### 2. cal Psi and Phi using gemma2 #######################################
###########################################################################
#!!! NOTE: we need to update tthis function using jax or C++/Rcpp with TG

# we cal Psi and Phi based on the null model (i.e., without X)

# Psi = matrix(c(2, -0.5, -0.5, 0.7), nrow = 2)
# Phi = matrix(c(0.35, 0.38,0.38, 0.76), nrow = 2)
res <- gemma2::MphEM(
  eval = K_eigen$values,
  X = t(matrix(rep(1, n_indi))) %*% K_eigen$vectors,
  Y = t(Y) %*% K_eigen$vectors,
  V_g = diag(ncol(Y)),
  V_e = diag(ncol(Y))
)[[1]]
res[c(2:3)]

Psi <- res$Vg
Phi <- res$Ve

Psi_eigen <- eigh(Psi)
Psi_Kva <- Psi_eigen$values
Psi_Kve <- Psi_eigen$vectors

Phi_eigen <- eigh(Phi)
Phi_Kva <- Phi_eigen$values
Phi_Kve <- Phi_eigen$vectors

Psi_Kva[Psi_Kva == 0] = 1e-13
Phi_Kva[Phi_Kva == 0] = 1e-13

###########################################################################
### 3. Data rotation ######################################################
###########################################################################

# Diagonalizing two matrices (Furlotte et al, Genetics 2015)
R <- Psi_Kve %*% sqrt(solve(diag(Psi_Kva))) # Now R %*% t(R) is equal to Psi^{-1}
RR <- (t(R) %*% Phi) %*% R
eigen_values <- eigh(RR)$values
D <- eigen_values
rightTransform <- t(eigh(RR)$vectors) %*% t(R)


# P: the variance of the jth phenotype befoe rotation (i.g., after transformation)
# P <- c(Kva + D[1], Kva + D[2], Kva + D[3])
P <- lapply(seq_len(length(D)), function(d){
  Kva + D[d]
}) %>% unlist

# The transformed vector of Y -> Yt / Y_T = vec(YM') ------->>>>>
Yt <- as.vector(Ystar %*% t(rightTransform)) %>% as.matrix()
#Yt <- as.vector(Ystar %*% rightTransform) %>% as.matrix()

# The transformed vector of Y = Yt
Xt <- kronecker(rightTransform, Xstar)


# rotation function using chol (Joo et al, Genetics 2016)
chol_solve <- function(K) {
  a = eigen(K)$vectors
  b = eigen(K)$values
  b[b < 1e-13] = 1e-13
  b = 1 / sqrt(b)
  return(a %*% diag(b) %*% t(a))
}

rotate <- function(Y, sigma) {
  U <- chol_solve(sigma)
  tU <-t(U)
  UY = tU%*%Y
  return(UY)
}

# data rotation 
new_x <- rotate(Xt, diag(P))
new_y <- rotate(Yt, diag(P))

###########################################################################
### 4-1. resutls single SNPs ##############################################
###########################################################################

# original results    
t1 <- solve(t(Xt) %*% solve(diag(P)) %*% Xt)
t2 <- t(Xt) %*% solve(diag(P)) %*% Yt
result <- t1 %*% t2
result
#           [,1]
# [1,] -0.689063
# [2,] -1.432046


single_anlaysis <- function(new_y, new_x){
  n_snp <- ncol(new_x) / 2

  res <- lapply(seq_len(n_snp), function(xx){
    tmpres <- summary(lm(new_y ~ new_x[, c(xx, (xx + n_snp))] - 1))$coefficients
    n_pheno <- nrow(tmpres)
    
    tmpres2 <- lapply(seq_len(n_pheno), function(pp) {
      sel_tmpres <- tmpres[pp, ]
      sel_tmpres <- t(data.frame(sel_tmpres))
      colnames(sel_tmpres) <- paste0(c("beta", "betastderr", "t", "p"), "_", pp)
      sel_tmpres %>% as.data.table
    })

    do.call(cbind.data.frame, tmpres2)
  }) %>% rbindlist
  return(res)
}

t1 <- solve(t(new_x) %*% new_x)
t2 <- t(new_x) %*% new_y
result <- t1 %*% t2
result

cbind(geno[c(sel_idx), c(1:3)], single_anlaysis(new_y, new_x))
#           X1 X2 X3     beta_1 betastderr_1       t_1          p_1     beta_2 betastderr_2       t_2
# 1  rs8275764  A  G -0.6890630    0.2852926 -2.415286 1.663145e-02 -1.4320457    0.3107636 -4.608152
#            p_2
# 1 7.268280e-06


cbind(geno[c(sel_idx), c(1:3)], single_anlaysis(new_y, new_x))
#           X1 X2 X3     beta_1 betastderr_1       t_1          p_1     beta_2 betastderr_2       t_2
# 1  rs8275764  A  G -0.6890630    0.2852926 -2.415286 1.663145e-02 -1.4320457    0.3107636 -4.608152
# 2  rs6212654  A  G -0.7908308    0.1889947 -4.184407 4.296073e-05  0.2302347    0.1806420  1.274536
# 3 rs13477740  G  A -0.6562880    0.2505862 -2.619011 9.501074e-03 -1.1481793    0.2659750 -4.316870
#            p_2
# 1 7.268280e-06
# 2 2.039674e-01
# 3 2.499081e-05


###########################################################################
### 4-2. resutls multiple snps ############################################
###########################################################################

gemma <- fread("GEMMA/output/mouse100.assoc.txt")
gemma %>% arrange(p_wald)
  #      chr            rs ps n_miss allele1 allele0    af        beta_1       beta_2  Vbeta_1_1   Vbeta_1_2
  #   1:  -9     rs8275764 -9      0       A       G 0.055 -6.683062e-01 -1.240216000 0.08150291 0.017920470
  #   2:  -9     rs8275844 -9      0       G       C 0.055 -6.683062e-01 -1.240216000 0.08150291 0.017920470
  #   3:  -9     rs6212654 -9      0       A       G 0.225 -7.578945e-01  0.198888000 0.03110525 0.007200907
  #   4:  -9    rs13477740 -9      0       G       A 0.070 -6.336982e-01 -0.993420600 0.06142864 0.011513500
  #   5:  -9 gnf04.060.296 -9      0       C       A 0.065 -7.059256e-01 -1.047843000 0.06685996 0.015020600
  #       Vbeta_2_2       p_wald
  #   1: 0.06438583 3.282434e-06
  #   2: 0.06438583 3.282434e-06
  #   3: 0.02684931 6.049233e-06
  #   4: 0.04851080 9.923329e-06
  #   5: 0.05717564 1.548445e-05
  #  ---                        