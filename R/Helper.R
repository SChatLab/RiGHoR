#### Per-Gene Modeling Function ####

perGene.Mckean <- function(expr, formula, regData){

  ### Fitting RLM ###

  regData <- data.frame(expr, regData)
  fit <- tryCatch({
    fit <- suppressWarnings(rfit(formula, data = regData))
  }, error=function(err){
    fit <- NA
    return(fit)
  })

  ### Collecting Outputs ###

  if(is.list(fit)){
    resDF <- as.data.frame(summary(fit)$coefficients)
    if(rownames(resDF)[1] == '1'){
      rownames(resDF) <- c('(Intercept)', expVar)
    }
    Beta <- resDF[grep(expVar, rownames(resDF)), 'Estimate']
    Intercept <- resDF['(Intercept)', 'Estimate']
    log2FC <- Beta
    se <- resDF[grep(expVar, rownames(resDF)), 'Std. Error']
    U.CI <- log2FC + qnorm(.025) * se
    L.CI <- log2FC - qnorm(.025) * se
    pval <- resDF[grep(expVar, rownames(resDF)), 'p.value']
  }else{
    log2FC <- NA
    se <- NA
    U.CI <- NA
    L.CI <- NA
    pval <- NA
  }

  ### Summarizing Outputs ###

  est.df <- data.frame(log2FC = log2FC,
                       SE = se,
                       L.CI = L.CI,
                       U.CI = U.CI,
                       Pval = pval)
  return(est.df)
}

#### Normalizing using TMM (edgeR) ####

tmm_norm <- function(features, metadata){
  norm.y <- DGEList(features)
  norm.y <- edgeR::calcNormFactors(norm.y, method = "TMM")
  norm.y <- as.data.frame(edgeR::cpm(norm.y, log = FALSE))
  return(norm.y)
}

#### Normalizing using Quantile ####

quan_norm <- function(features,metadata){
  norm.y <- normalize.quantiles(as.matrix(features))
  norm.y <- data.frame(norm.y)
  names(norm.y) <- names(features)
  rownames(norm.y) <- rownames(features)
  return(norm.y)
}

#### Normalizing using Upper Quartile ####

uqrt_norm <- function(features,metadata){
  quant.exp <- apply(as.matrix(features), 2, function(x){quantile(x[x > 0], 0.75)})
  norm.y <- data.frame(t(t(as.matrix(features))/quant.exp))
  return(norm.y)
}

#### Normalizing using Geometric Means (DESeq2) ####

rle_norm <- function(features, metadata, coVars = NULL, expVar = 'Exposure'){
  if(is.null(coVars)){
    metadata <- metadata[, c(expVar), drop = FALSE]
  }else{
    metadata <- metadata[, c(expVar, coVars)]
  }
  formula <- as.formula(paste('~', paste(colnames(metadata), collapse = "+"), sep = ''))
  x <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(features),
                                                       colData = metadata,
                                                       design = formula))
  gm_mean <- function(x, na.rm = TRUE){
    exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
  }
  geoMeans <- apply(DESeq2::counts(x), 1, gm_mean)
  s <- DESeq2::estimateSizeFactors(x,geoMeans = geoMeans)
  s <- s$sizeFactor
  norm.y <- data.frame(t(apply(features, 1, function(x)x/s)))
  return(norm.y)
}
