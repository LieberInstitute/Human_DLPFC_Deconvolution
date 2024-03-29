
run_DE <- function(rse, model, run_voom = TRUE, save_eBayes = FALSE, coef, plot_name = NULL){
  
  print_plots <- !is.null(plot_name)
  
  if(print_plots) pdf(plot_name)
  
  ## limma
  eBayes_out = get_eBayes(rse = rse, model = model, run_voom = run_voom, print_plots)
  
  #### because of more than 1 component, computing F statistics instead of t=-statistics
  # message("Calc top tables - coefficents:", paste(coef, collpase = " "))
  if(all(is.character(coef))){
    message("Character coef in eBayes: ", all(coef %in% colnames(eBayes_out$coefficients)))
  }else{
    message("Index coef: ", colnames(eBayes_out$design)[[coef]])
    coef <- colnames(eBayes_out$design)[[coef]]
  }
  
  
  topTable_out = topTable(eBayes_out, coef=coef, number=Inf , sort.by = "none")
  
  if(print_plots){
    message("plotting: ", plot_name)
    try(plotMA(eBayes_out, coef = coef))
    try(hist(topTable_out$P.Value))
    dev.off()
  }
  
  if(length(coef == 1)){
    message("1 coef, Return t-stats")
    topTable_out[,paste0("q_", coef)] <- p.adjust(topTable_out$P.Value, 'fdr')
    
  } else {
    message("Multiple coef, Return F-stat")
    ## significance levels EXTRACT INDIVIDUAL COMPARISON P-VALUES THAT ARE NOT IN TOP TABLE
    pvalMat = as.matrix(eBayes_out$p.value)[,coef]
    qvalMat = pvalMat
    qvalMat[,1:2] = p.adjust(pvalMat[,1:2],method="fdr") ## TODO fix 1:2 to stable names
    colnames(pvalMat) = paste0("P_",colnames(pvalMat))
    colnames(qvalMat) = paste0("q_",colnames(qvalMat))
    
    topTable_out = cbind(topTable_out,cbind(pvalMat, qvalMat))
  }
  
  ## print summary
  n_pVal <- report_top_pVal(topTable_out = topTable_out, cols = paste0("q_", coef))
  print(n_pVal)
  
  if(save_eBayes){
    output <- list(topTable = topTable_out,
                   eBayes = eBayes_out)
    return(output)
  } else {
    return(topTable_out)
  }
  
}

get_eBayes <- function(rse, model, run_voom = TRUE, print_plots = FALSE){
  
  if(run_voom){
    message(Sys.time(), " - Calc Norm Factors")
    dge_out = DGEList(counts = assays(rse)$counts,
                      genes = rowData(rse))
    dge_norm = calcNormFactors(dge_out)
    voom_out = voom(dge_norm, model, plot=print_plots)
    
    message(Sys.time(), " - Limma")
    fit = lmFit(voom_out)
  } else {
    message("Using tpm")
    log_tmp = log2(assays(rse)$tpm + 1)
    
    message(Sys.time(), " - Limma")
    fit = lmFit(log_tmp, model)
  }
  
  ## limma
  message(Sys.time(), " - eBayes")
  eBayes_out = eBayes(fit)
  return(eBayes_out)
}

report_top_pVal <- function(topTable_out, cols){
  
  tt_values <- topTable_out[,cols, drop = FALSE]
  cutoffs <- list(0.05, 0.01)
  names(cutoffs) <- paste("< ", cutoffs)
  n_sig <- map(cutoffs, ~colSums(tt_values < .x))
  sig_table <- t(do.call("rbind", n_sig))
  return(sig_table)
}
