#- Nikos: you can use a code like above! where measurements is basically you X. And this will give you 2 matrices.
#  One with rows the genesets and columns the sample that has Normalized Enrichment Scores (NES) and one with the same 
#  dimesnsions and p-values (for the later we dont care at this point so you can avoid using it if you want and for that
#   reason you may use just 100 permutations)

# code for fgsea from nikos. get a score for every sample in the dataset
print("running fgsea for enrichment space")
  genesets_list <-apply(measurements,MARGIN = 2,fgsea,pathways = all_genesets,
                  minSize=10,
                  maxSize=500,
                  nperm = n_permutations)
  print("fgsea finished, preparing outputs")
  
  # Prepare output
  
  NES <- genesets_list[[1]]$NES
  if (pval_adjustment==T){
    pval <- genesets_list[[1]]$padj
  }else{
    pval <- genesets_list[[1]]$pval
  }
  if (length(genesets_list)>1){
    for (i in 2:length(genesets_list)) {
      
      NES <- cbind(NES,genesets_list[[i]]$NES)
      if (pval_adjustment==T){
        pval <- cbind(pval,genesets_list[[i]]$padj)
      }else{
        pval <- cbind(pval,genesets_list[[i]]$pval)
      }
    }
  }else{
    NES <- as.matrix(NES)
    pval <- as.matrix(pval)
  }
  
  colnames(NES) <- names(genesets_list)
  rownames(NES) <- genesets_list[[1]]$pathway
  colnames(pval) <- names(genesets_list)
  rownames(pval) <- genesets_list[[1]]$pathway
