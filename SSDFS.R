rm(list=ls())
options(stringsAsFactors = F)
setwd("./")
getwd()

library(cluster)
library(tidyverse)

rnaseq_processing <- function(data, log2 = T) {
  
  pos <-
    which(apply(data, 1, function(x) {
      length(which(x == 0))
    }) < 0.94 * ncol(data))
  data <- data[pos, ]
  
  pos <-
    which(apply(data, 1, function(x) {
      length(which(x == 0))
    }) > 0.06 * ncol(data))
  data <- data[pos, ]
  
  if (log2) {
    data <-  log2(data + 1)
  }
  return(data)
}

evaluation <- function(truelabel, predlabel) {
  if (length(truelabel) != length(predlabel))
    stop("truelabel and predlabel must have the same length")
  total = length(truelabel)
  x_ids = unique(truelabel)
  y_ids = unique(predlabel)
  MI = 0.0
  for (idx in x_ids) {
    for (idy in y_ids) {
      idxOccur = which(truelabel == idx)
      idyOccur = which(predlabel == idy)
      idxyOccur = intersect(idxOccur, idyOccur)
      if (length(idxyOccur) > 0) {
        MI = MI + (length(idxyOccur) / total) * log2((length(idxyOccur) * total) /
                                                       (length(idxOccur) * length(idyOccur)))
      }
    }
  }
  Hx = 0
  for (idx in x_ids) {
    idxOccurCount = length(which(truelabel == idx))
    
    Hx = Hx - (idxOccurCount / total) * log2(idxOccurCount / total)
  }
  Hy = 0
  for (idy in y_ids) {
    idyOccurCount = length(which(predlabel == idy))
    Hy = Hy - (idyOccurCount / total) * log2(idyOccurCount / total)
  }
  nmi = 2 * MI / (Hx + Hy)
  tab = table(truelabel, predlabel)
  conv_df = as.data.frame.matrix(tab)
  n <- sum(tab)
  ni <- apply(tab, 1, sum)
  nj <- apply(tab, 2, sum)
  n2 <- choose(n, 2)
  nis2 <- sum(choose(ni[ni > 1], 2))
  njs2 <- sum(choose(nj[nj > 1], 2))
  ri = 1 + (sum(tab ^ 2) - (sum(ni ^ 2) + sum(nj ^ 2)) / 2) / n2
  ari = c(sum(choose(tab[tab > 1], 2)) - (nis2 * njs2) / n2) / ((nis2 + njs2) /
                                                                  2 - (nis2 * njs2) / n2)
  out = c(nmi,ari)
  names(out)=c("NMI","ARI")
  return(out)
}

HC_clustering <- function(dataset, K) {
  distance <- as.dist(1 - cor(dataset))
  h_tree <- hclust(distance, method = "complete")
  cluster_result <- cutree(h_tree, k = K)
  return(cluster_result)
}


load("./dataset/Goolam.Rdata")
label <- label
table(label)
k_n <- length(unique(label))
data_set <- rnaseq_processing(data)
data_set <- t(data_set)

X <- data_set
XX <- t(X) %*% X
XXX <- XX %*% XX
n <- nrow(X)
d <- ncol(X)

NITER <- 50
k <- 500
rho <- 1
rho1 <- 1
set.seed(666)
W <- matrix(runif(d * k), nrow = d, ncol = k)
obj <- numeric(NITER)
feature_set <- vector("list", length = NITER)

for (iter in 1:NITER) {
  cat("<<<Number of iterations: ", iter, "\n")
  numW <- XXX %*% W + rho * W
  A <- XX %*% W
  denW <- A %*% t(W) %*% A + rho * matrix(1, nrow = d, ncol = d) %*% W + rho1 * XX %*% W %*% matrix(1, nrow = k, ncol = k)
  fracW <- numW / denW
  fracW <- fracW^(1/4)
  W <- W * fracW
  XW <- X %*% W
  WW <- W %*% t(W)
  obj[iter] <- 0.5 * norm(X %*% t(X) - XW %*% t(XW), 'F')^2 + rho * (sum(diag(matrix(1, nrow = d, ncol = d) %*% WW)) - sum(diag(WW))) + rho1 * matrix(1, nrow = 1, ncol = k) %*% t(XW) %*% XW %*% matrix(1, nrow = k, ncol = 1)           
  cat("<<<objective function: ",obj[iter],"\n")
  score <- rowSums(W * W)
  index <- order(score, decreasing = TRUE)
  feature_set[[iter]] <- index[1:k]
  cat("####################################","\n")
}

feature_set
