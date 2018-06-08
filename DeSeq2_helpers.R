#Helper functions

#This one simulates a single DESeq experiment and analyses its statistics
deSeqRun <- function(percent_with_effect, effect_size, lfcThreshold, num_replicates, num_genes) {
  sim_data <- simulateDeSeq2(num_genes, num_replicates, percent_with_effect, effect_size)
  dds <- DESeqDataSetFromMatrix(countData = sim_data$observed$counts,
                                colData = sim_data$observed$sample_info,
                                design= ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("condition","1","2"), lfcThreshold = lfcThreshold)
  resultsNames(dds)
  res_shrink <- lfcShrink(dds, coef=2)
  
  sum_func <- function(x) { sum(x, na.rm  = TRUE)}
  res_df <- res %>% as.data.frame() %>%
    mutate(
      log2FoldChange_shrunk = -res_shrink$log2FoldChange, #Note the minus with res_shrink - for some reason lfcShrink inverts the sign :-O
      true_effect = sim_data$true$coefficients[,1], truly_changed = true_effect != 0, significant = padj < 0.05
    ) %>%
    mutate(TP = truly_changed & significant, TN = !truly_changed & !significant, 
           FP = !truly_changed & significant, FN = truly_changed & !significant,
           S_error = if_else(significant, true_effect * log2FoldChange < 0, FALSE)
    )
  
  stats <- res_df %>%
    summarise(TP_ = sum_func(TP), 
              TN_ = sum_func(TN), 
              FP_ = sum_func(FP), 
              FN_ = sum_func(FN), 
              S_error_ = sum_func(S_error), 
              mean_true_eff = mean(if_else(truly_changed & significant, abs(log2FoldChange), as.numeric(NA)), na.rm = TRUE),
              mean_true_eff_shrunk = mean(if_else(truly_changed & significant, abs(log2FoldChange_shrunk), as.numeric(NA)), na.rm = TRUE),
              mean_false_eff = mean(if_else(!truly_changed & significant, abs(log2FoldChange), as.numeric(NA)), na.rm = TRUE),
              median_p_positive = median(padj[significant], na.rm = TRUE),
              median_p_negative = median(padj[!significant], na.rm = TRUE),
              median_p_fn = median(padj[FN], na.rm = TRUE)
    ) %>%
    mutate(num_replicates = num_replicates, percent_with_effect = percent_with_effect, effect_size = effect_size, lfcThreshold = lfcThreshold)
  stats
  #data.frame(true = sim_data$true$coefficients[,2], estimated = res$baseMean) %>% filter(true != 0 & estimated != 0) %>% mutate(estimated = log(estimated)) %>% ggplot(aes(x=true, y = estimated)) + geom_point()
  
}

#Runs multiple DESeq runs and gathers the statistics
deSeqTest <- function(percent_with_effect, effect_size, lfcThreshold, num_replicates = 3, num_genes = 1000, num_simulations = 100) {
  apply(array(1:num_simulations), 1, FUN = function(x) {deSeqRun(percent_with_effect, effect_size,  lfcThreshold, num_replicates,num_genes)}) %>% do.call(rbind,.)  
}

deSeqMultiTest <- function(inputs, num_simulations = 100) {
  num_simulations = 100
  results_base_list <- list()
  for(i in 1:nrow(inputs)) {
    results_base_list[[i]] <- deSeqTest(0.2, effect_size = inputs$effect_size[i], lfcThreshold = inputs$lfcThreshold[i], num_simulations = num_simulations, num_replicates = inputs$num_replicates[i])
  }
  
  do.call(rbind, results_base_list)
}


#This function simulates one replicated experiment (basically simulates one larger experiemtn and splits it into two)
deseq_replication <- function(percent_with_effect, effect_size, lfcThreshold, num_replicates = 3, num_genes = 1000) {
  
  sim_data_both <- simulateDeSeq2(num_genes, num_replicates * 2, percent_with_effect, effect_size)
  counts1 <- seq(1,num_replicates * 4, by = 2)
  counts2 <- seq(2,num_replicates * 4, by = 2)
  dds1 <- DESeqDataSetFromMatrix(countData = sim_data_both$observed$counts[,counts1],
                                 colData = sim_data_both$observed$sample_info %>% filter(id %in% counts1),
                                 design= ~ condition)
  dds1 <- DESeq(dds1)
  res1 <- results(dds1, contrast=c("condition","1","2"),lfcThreshold = lfcThreshold)
  
  dds2 <- DESeqDataSetFromMatrix(countData = sim_data_both$observed$counts[,counts2],
                                 colData = sim_data_both$observed$sample_info %>% filter(id %in% counts2),
                                 design= ~ condition)
  dds2 <- DESeq(dds2)
  res2 <- results(dds2, contrast=c("condition","1","2"),lfcThreshold = lfcThreshold)
  
  res_df <- res1 %>% as.data.frame() %>% mutate(padj2 = res2$padj, lfc2 = res2$log2FoldChange) 
  
  sum_func <- function(x) { sum(x, na.rm  = TRUE)}
  
  res_df %>% mutate(significant = padj < 0.05) %>% summarise(
    significant_ = sum_func(significant), 
    replicated =  sum_func(padj < 0.05 & padj2 < 0.05),
    effect_dev = mean(lfc2 - log2FoldChange, na.rm = TRUE),
    significant_effect_dev = mean(lfc2[significant] - log2FoldChange[significant], na.rm = TRUE),
    smaller_eff = sum(abs(lfc2) < abs(log2FoldChange), na.rm = TRUE),
    smaller_eff_significant = sum(abs(lfc2[significant]) < abs(log2FoldChange[significant]), na.rm = TRUE)
  ) %>% 
    mutate(effect_size = effect_size, percent_with_effect = percent_with_effect, 
           lfcThreshold = lfcThreshold, num_replicates = num_replicates)
}

deseq_replication_multi <- function(inputs, num_simulations = 100) {
  results_list_repl <- list()
  j = 1
  
  
  for(i in 1:nrow(inputs)) {
    for(sim in 1:num_simulations) {
      results_list_repl[[j]] <- deseq_replication(0.2, effect_size = inputs$effect_size[i], lfcThreshold = inputs$lfcThreshold[i])
      j <- j + 1
    }
  }
  
  do.call(rbind, results_list_repl)
}