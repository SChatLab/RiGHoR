fit.RoNseq <- function(features, 
                       metadata,
                       expVar = 'Exposure', 
                       coVars = NULL,
                       norm.method = 'tmm',
                       parallel = FALSE,
                       ncores = 4){
  
  #### Creating Regression Pre-requisites ####
  
  start.time <- Sys.time()
  mem <- peakRAM::peakRAM({
    if(is.null(coVars)){
      regData <- metadata[, c(expVar), drop = FALSE]
    }else{
      regData <- metadata[, c(expVar, coVars)]
    }
    regData[sapply(regData, is.character)] <- lapply(regData[sapply(regData, is.character)], as.factor)
    formula <- as.formula(paste("expr ~ ", paste(colnames(regData), collapse = "+")))
    
    #### Normalizing Expression counts ####
    
    if(norm.method == 'tmm'){
      norm.y = log(suppressMessages(tmm_norm(features, metadata)) + 1)
    }else if(norm.method == 'rle'){
      norm.y = suppressMessages(rle_norm(features, metadata))
    }else if (norm.method == 'cpm'){
      norm.y <- cpm_norm(features = as.matrix(features), metadata, 
                         const.mult = 1e+06, prior.count = 1, log = FALSE)
    }else if (norm.method == 'uquart'){
      norm.y = uqrt_norm(features, metadata)
    }else if (norm.method == 'quant'){
      norm.y = quan_norm(features, metadata)
    }
    
    #### Apply model ####
    if(parallel){
      cl <- parallel::makeCluster(ncores)
      doSNOW::registerDoSNOW(cl)
      packages <- c('Rfit')
      exports <- c('perGene.Mckean', 'expVar', 'coVars')
      pb <- txtProgressBar(max = nrow(norm.y), style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      res <- foreach(j = 1:nrow(norm.y), .combine = rbind, 
                     .packages = packages, .options.snow = opts,
                     .export = exports) %dopar% {
                       expr <- as.numeric(norm.y[j, ])
                       tmpfit <- perGene.Mckean(expr = expr, 
                                                formula = formula, 
                                                regData = regData)
                       tmpfit <- data.frame(Gene = j, tmpfit)
                       return(tmpfit)
                     }
      close(pb)
      stopCluster(cl)
      rownames(res) <- rownames(norm.y)
    }else{
      res <- apply(norm.y, 1, function(expr) perGene.Mckean(expr = expr, 
                                                            formula = formula, 
                                                            regData = regData))
      res <- do.call('rbind', res)
    }
  })

  #### Summarizing Output ####

  Genes <- data.frame(Genes = rownames(features))
  res$Genes <- rownames(res)
  res$adjPval <- p.adjust(res$Pval, method = 'BH')
  output <- merge(Genes, res, id = 'Genes', all = TRUE)
  output <- output[order(output$adjPval), ]
  peak.memory.gb <- mem$Peak_RAM_Used_MiB * (2^20) / (10^9)
  stop.time <- Sys.time()
  Time.min = round(difftime(stop.time, start.time, units = 'mins')[[1]], 3)
  finres = list(Method = 'RoNseq',
                res = output, 
                Time.min = Time.min,
                peak.memory.gb = peak.memory.gb)
  return(finres)
}
