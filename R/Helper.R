#### Per-Gene Modeling Function ####

perGene.Mckean <- function(expr, formula, regData){
  
  ### Fitting RLM ###
  
  regData <- data.frame(expr, regData)
  
  safe_rfit <- function(scores_fun){
    tryCatch(
      suppressWarnings(rfit(formula, data = regData, scores = scores_fun)),
      error = function(e) NA
    )
  }
  
  # simple moment-based skewness and kurtosis
  skewness_m <- function(x){
    x <- x[is.finite(x)]
    if (length(x) < 5) return(NA_real_)
    m <- mean(x); s <- sd(x)
    if (s == 0 || is.na(s)) return(NA_real_)
    mean((x - m)^3) / (s^3)
  }
  
  kurtosis_m <- function(x){
    x <- x[is.finite(x)]
    if (length(x) < 5) return(NA_real_)
    m <- mean(x); s <- sd(x)
    if (s == 0 || is.na(s)) return(NA_real_)
    mean((x - m)^4) / (s^4)
  }
  
  pick_scores <- function(sk, kt){
    # If diagnostics are missing, fall back to wscores
    if (is.na(sk) || is.na(kt)) return(wscores)
    
    if (sk > 0.7) return(bentscores1)   # right-skew
    if (sk < -0.7) return(bentscores3)  # left-skew
    
    # tails
    if (kt < 3.0) return(bentscores2)   # light tails
    if (kt > 4.5) return(bentscores4)   # heavy tails
    
    # close-ish to normal
    if (kt >= 3.0 && kt <= 4.5) return(nscores)
    
    wscores
  }
  
  # baseline fit
  fit0 <- safe_rfit(wscores)
  
  # compute residual shape and pick a score
  if (is.list(fit0)) {
    rs <- tryCatch(rstudent(fit0), error = function(e) NA)
    sk <- skewness_m(rs)
    kt <- kurtosis_m(rs)
    
    scores <- pick_scores(sk, kt)
    fit <- safe_rfit(scores)
    
    diag_skew <- sk
    diag_kurt <- kt
  } else {
    fit <- NA
    scores <- NA
    diag_skew <- NA
    diag_kurt <- NA
  }
  
  ### Collecting Outputs ###
  
  if (is.list(fit)) {
    resDF <- as.data.frame(summary(fit)$coefficients)
    
    if (rownames(resDF)[1] == "1") {
      rownames(resDF) <- c("(Intercept)", expVar)
    }
    
    Beta <- resDF[grep(expVar, rownames(resDF)), "Estimate"]
    Intercept <- resDF["(Intercept)", "Estimate"]
    
    log2FC <- Beta
    se <- resDF[grep(expVar, rownames(resDF)), "Std. Error"]
    
    z <- qnorm(0.975)
    L.CI <- log2FC - z * se
    U.CI <- log2FC + z * se
    
    pval <- resDF[grep(expVar, rownames(resDF)), "p.value"]
  } else {
    log2FC <- NA
    se <- NA
    U.CI <- NA
    L.CI <- NA
    pval <- NA
    Beta <- NA
    Intercept <- NA
  }
  
  ### Summarizing Outputs ###
  
  est.df <- data.frame(log2FC = log2FC, 
                       SE = se,
                       L.CI = L.CI,
                       U.CI = U.CI,
                       Pval = pval)
  return(est.df)
}
