
mod.sumStat <- function(mm = NULL, data_list = NULL, poiss = FALSE){
  
  for(i in 1:length(data_list)){
    assign(names(data_list)[i], data_list[[i]])
  }
  
  lik_mat.det <- matrix(0, ncol = T*K*M, nrow = nrow(mm))
  lik_mat.eco <- matrix(0, ncol = 2*(T-1)+1, nrow = nrow(mm))
  
  cat(paste0("\n--------------------------------------------------| ", nrow(mm), "\n"))
  for(o in 1:nrow(mm)){
    if(o%%floor(nrow(mm)/50) == 0 ){cat("+")}
    a <- mm[o,grep("a", colnames(mm))]
    b <- mm[o,grep("b", colnames(mm))]
    d <- mm[o,grep("d", colnames(mm))]
    lam <- mm[o,"lam"]
    
    a_RE_stat <- mm[o,grep("a_RE_stat", colnames(mm))]
    a <- a[1:Ncov_rho]
    
    n <- mm[o,grep("n", colnames(mm))]
    G <- mm[o,grep("G", colnames(mm))]
    S <- mm[o,grep("S", colnames(mm))]
    
    if(poiss){
      rho <- foreach(m = 1:M, .combine = rbind)%do%{
        exp(a_RE_stat[m] + a %*% t(rho_covs))}
    }
    if(!poiss){
      rho <- foreach(m = 1:M, .combine = rbind)%do%{
        invlogit(a_RE_stat[m] + a %*% t(rho_covs))}
    }
    
    gam <- exp(b %*% t(gam_covs))
    tau <- invlogit(d %*% t(tau_covs))
  
  lik.det <- array(0, dim = c(M, T, K))
  lik.eco <- rep(0, 2*(T-1)+1)
  lik.eco[1] <- dpois(n[1], lam)
  for(t in 1:T){
    if(t != T){
      lik.eco[1+t] <- dpois(G[t], gam[t])
      lik.eco[T+t] <- dbinom(S[t], n[t],tau[t])
      }
    for(m in 1:M){
      p <- 1 - (1 - rho[m,t])**n[t]
      for(k in 1:K){
        if(poiss){lik.det[m,t,k] <- dpois(y[m,t,k], n[t]*rho[m,t])}
        if(!poiss){lik.det[m,t, k] <- c(1-p, p)[y[m,t,k]]}
      }
    }
  }
  lik_mat.det[o,] <- as.numeric(lik.det)
  lik_mat.eco[o,] <- as.numeric(lik.eco)
  }
  CPO.det <- - sum(log(nrow(mm) / apply(1/lik_mat.det,2,function(x){sum(x)})), na.rm = T)
  CPO.eco <- - sum(log(nrow(mm) / apply(1/lik_mat.eco,2,function(x){sum(x)})), na.rm = T)

  WAIC.det <- -2 * sum(log(apply(lik_mat.det,2,function(x){mean(x)})), na.rm = T) +
    sum(apply(lik_mat.det,2,function(x){var(x)}), na.rm = T)
  WAIC.eco <- -2 * sum(log(apply(lik_mat.eco,2,function(x){mean(x)})), na.rm = T) +
    sum(apply(lik_mat.eco,2,function(x){var(x)}), na.rm = T)
  
  list(CPO.det = CPO.det, CPO.eco = CPO.eco, CPO = CPO.det + CPO.eco,
       WAIC.det = WAIC.det, WAIC.eco = WAIC.eco, WAIC = WAIC.det + WAIC.eco)
}