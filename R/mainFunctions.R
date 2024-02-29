#' hiernet.lasso.binary
#' This function adapts the hierNet::hierNet.path() function such that it has the correct
#' input format for the stabs::stabsel() function (similar to
#' stabs::glmnet.lasso(). As hierNet::hierNet.path() initially
#' was not designed for binary input features the internal scaling is inappropriate 
#' and zz.binary needs to be computed and passed to hierNet::hierNet.path().
#' @param x binary input matrix of the dimension n x p (samples x proteins)
#' @param y outcome vector of length n
#' @param q number of (unique) selected variables (or groups of variables depending on the model) that are selected on each subsample.
#' @param l number of main effects/linear features -> stabsel() requires all effects
#'         (main + interactions) and hierNet.path() only requires linear features as x input


library(hierNet)
library(stabs)
library(parallel)

# read data

# L <- readRDS("L-library.rds")
# P <- readRDS("P-proteinsFandR.rds")



hiernet.lasso.binary <- function(x, y, q, l = ncol(L_p)#, strong_ = TRUE
){
  
  if (!requireNamespace("hierNet", quietly = TRUE))
    stop("Package ", sQuote("hierNet"), " needed but not available")
  
  
  x.lin <- x[, 1:l]
  p_init <- ncol(x.lin)
  ## there might be zero columns (modifications) after removing observations
  ## (subsampling) which compute.interactions can't handle:
  ## Remove zero columns
  col_sub = apply(x.lin, 2, function(col) all(col == 0 ))
  x.lin = x.lin[, !col_sub]
  
  
  nlam = 20
  
  
  
  
  ## with passing zz we avoid that zz gets computed based on the scaled x
  ## which is not wanted for binary features
  zz_binary <- compute.interactions.c(x.lin, diagonal = FALSE)
  fit <- hierNet::hierNet.path(x.lin, y, minlam = .001, maxlam = 20, 
                               nlam = nlam, diagonal = FALSE, strong = TRUE,#strong_,
                               stand.int = FALSE, stand.main = TRUE,
                               zz = zz_binary
  )
  
  
  
  
  p <- ncol(x.lin)
  coefmatrix_lam <- matrix(nrow = p * (p + 1)/2, ncol = nlam)
  
  linear_coef <- fit$bp - fit$bn ## linear coefficients
  
  for(lam in 1:length(fit$lamlist)){
    
    theta <- (fit$th[, , lam] + t(fit$th[, , lam]))/2 ## symmetrize theta
    interact_coef <- theta[lower.tri(theta)]
    
    coefmatrix_lam[, lam] <- c(linear_coef[, lam], interact_coef) ## combine linear effects and interactions
    
  }
  
  ## interaction feature names after potentially removing linear features
  A <- colnames(x.lin)
  n <- 0
  names_int <- c()
  for(i in 1:(ncol(x.lin) - 1)){
    for(j in (i + 1): ncol(x.lin)){
      n <- n + 1
      names_int[n] <- paste0(A[i], ":", A[j])
    }
  }
  
  rownames(coefmatrix_lam) <- c(colnames(x.lin), names_int)
  
  ## Empty matrix with all features (also the ones that where removed before)
  sequence <- matrix(nrow = p_init * (p_init + 1)/2, ncol = nlam)
  rownames(sequence) <- colnames(x)
  sequence[rownames(coefmatrix_lam), ] <- (coefmatrix_lam != 0)
  sequence[is.na(sequence)] <- 0
  
  ## select the lambda where number of selected features is <= q
  seq_q <- which(colSums(sequence) <= q)
  ret <- sequence[, seq_q[length(seq_q)]]
  
  ## return both
  return(list(selected = ret, path = sequence))
}

P_ <- P
L_ <- L



#' stabsel_hiernet
#' This function runs hiernet with stability selection
#
#' @param p p-th feature in P n x p (samples x proteins)
#' @param P n x p - dimensional outcome
#' @param L n x q experimental design





stabsel_hiernet <- function(p, P = P_, L = L_){
  
  
  p <- unlist(p)
  
  
  
  ## indices of sign inconsistent measurements (F vs. R)
  # sign_incon <- which(sign(as.vector(P[, p])) !=
  #                       sign(as.vector(P[, p + 1915])))
  
  sign_incon <- which(abs(as.vector(P[, p]) - as.vector(P[, p + 1915])) > 0.3)
  ## remove sign inconsistent observations
  if(length(sign_incon) != 0){
    PF_p <- P[-sign_incon, p] ## forward protein binding
    PR_p <- P[-sign_incon, p + 1915] ## reverse protein binding
    L_p <- L[-sign_incon, ] ## Library
  }
  else{
    PF_p <- P[, p]
    PR_p <- P[, p + 1915]
    L_p <- L
  }
  
  ## Remove zero columns
  col_sub = apply(L_p, 2, function(col) all(col == 0 ))
  L_p = L_p[, !col_sub]
  L_p_intactions <- cbind(L_p, hierNet::compute.interactions.c(L_p, diagonal = F)) ## add interactions
  
  
  hiernet.lasso.binary_ <- function(x, y, q, l = ncol(L_p)) { return(hiernet.lasso.binary(x, y, q, l = ncol(L_p))) }
  
  ## Fit model Forward experiment
  set.seed(3456)
  stab.hiernetF <- stabsel(x = L_p_intactions, y = PF_p,
                           fitfun = hiernet.lasso.binary_, cutoff = 0.6,
                           q = 12, sampling.type = "SS")
  ## Fit model Reverse experiment
  set.seed(3456)
  stab.hiernetR <- stabsel(x = L_p_intactions, y = PR_p,
                           fitfun = hiernet.lasso.binary_, cutoff = 0.6,
                           q = 12, sampling.type = "SS")
  
  return(list("F" = stab.hiernetF, "R" = stab.hiernetR))
}

# stabsel_hiernet_cluster_dist <- parallel::mclapply(as.list(1:1915),
#                                               function(p, P, L){return(stabsel_hiernet(p, P, L))},
#                                               P = P_, L = L_,
#                                               mc.cores = 8, mc.preschedule = FALSE)
# 
# saveRDS(stabsel_hiernet_cluster_dist, "/Users/mara.stadler/LRZ Sync+Share/PhD/MUDS_Projekt/temp/stabsel_hiernet_cluster__dist.rds")




run_workflow <- function(L, P, CPSS_F, CPSS_R, cutoff = .5){
  keep_or_not <- c()
  ## how many proteins have exaclty the same model in F and R:
  same_model <- c()
  n <- 0
  
  L_interactions <- cbind(L, hierNet::compute.interactions.c(L, diagonal = F)) ## add interactions
  
  pp <- ncol(P)/2
  fit_lm <- list()
  coef_matrix <- matrix(NA, nrow = ncol(L_interactions), ncol = pp)
  rownames(coef_matrix) <- colnames(L_interactions)
  
  
  
  ## binary matrix showing selected features after filtering step 2
  sel_matrix_FR <- matrix(0, nrow = ncol(L_interactions), ncol = pp)
  rownames(sel_matrix_FR) <- colnames(L_interactions)
  colnames(sel_matrix_FR) <- substr(colnames(P[, 1:pp]), 1,
                                    nchar(colnames(P[, 1:pp]))-2)
  
  for(p in 1:pp){
    if(class(CPSS_F[[p]]) != "empty" | class(CPSS_R[[p]]) != "empty" ){
      
      ## define set of selected features for both replicates
      selF_p <- which(CPSS_F[[p]]$max >= cutoff)
      selR_p <- which(CPSS_R[[p]]$max >= cutoff)
      
      ## only do refitting on mean over replicates if selected features overlap 
      ## or are a superset of each other
      if(length(selF_p) != 0 & length(selR_p) != 0){
        if(length(intersect(selR_p, selF_p)) > 0){
          if(all(selF_p %in% selR_p) | all(selR_p %in% selF_p)){
            
            if(all(selF_p %in% selR_p) & all(selR_p %in% selF_p)){
              same_model[p] <- 1}
            
            
            keep_or_not[p] <- 1
            n <- n + 1
            Pmean_p <- log((2^P[, p] + 2^P[, p + pp])/2, base = 2)
            Lint <- L_interactions[, intersect(selR_p, selF_p)]
            sel_matrix_FR[intersect(selR_p, selF_p), p] <- 1
            
            fit_lm[[p]] <- lm(Pmean_p ~ Lint)
            beta <- fit_lm[[p]]$coefficients[-1]
            names(beta) <- colnames(Lint)
            coef_matrix[names(beta), p] <- beta
            
            
          }
          else{keep_or_not[p] <- 0}
        }
      }
    }
    
    else{keep_or_not[p] <- 0}
    
  }
  colnames(coef_matrix) <- substr(colnames(P[, 1:1915]), 1,
                                  nchar(colnames(P[, 1:1915]))-2)
  ## set remaining (not selected) values to zero
  coef_matrix[is.na(coef_matrix)] <- 0
  return(list("coef_matrix" = coef_matrix, "keep_or_not" = keep_or_not, "refit_lm" = fit_lm))
}
