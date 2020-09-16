load("data/appdata.RData")
library(sjmisc)
library(dplyr)
library(FoldGO)


query <- function(varlist, smodel){
  out.ids <- list()
  model <- switch(smodel,
                  "Fisher's pvalue" = meta_comb@metaresult,
                  "Random Effects Modeling" = meta_degs_rem@metaresult)
  
  
  for(i in 1:length(varlist)){
    #getting the intersection of queries from the REM dataset
    if(!(is.null(varlist[i])) && smodel == "Random Effects Modeling"){
      name <- names(varlist)[i]
      out.var <- switch(name,"pval" = model["randomP"] <= varlist[i],
                        "logfc2" = model["randomSummary"] >= varlist[i],
                        "tRank" = model["rank"] <= varlist[i])
      
      out.ids[name] <- list(out.var)
    }
    
    
    
    
    else if(!(is.null(varlist[i])) && smodel == "Fisher's pvalue"){
      out.var <- switch(as.character(i), "pval" = (model["metaP"] <= pval),
                        "logfc2" = model["metafc"] >= fc)
      
      out.ids[i] <- list(out.var)
      return(intersect(out.ids[1][1], out.ids[2][1]))
      
    }
    
  }
  dfs <- lapply(out.ids, data.frame, stringsAsFactors = FALSE)
  binded.ids <- bind_cols(dfs)
  rownames(binded.ids) <- model$Symbol
  
  return(row.check(binded.ids))
}
row.check <- function(df){
  #rotate my data.frame bb
  df <- t(df)
  #output
  out <- c()
  #iterate and find the rows that are all true
  for(col in 1:length(df[1,])){
    if(all(df[,col])){
      out <- append(out, colnames(df)[col]) 
    }
  }
  return(out)
  
}

median_dups <- function(logfc2, dupList){
  medians <- list()
  for(i in 1:length(dupList)){
    fc <- median(logfc2[as.vector(dupList[[i]])])
    medians[names(dupList)[i]] <- list(fc)
  }
  return(medians)
}

make_dup_list <- function(vnames){
  out <- list()
  i <- 1
  q_list <- vnames
  while(length(q_list > 1)){
    
    query <- vnames[i]
    
    q <- which(vnames == query)
    out[query] <- list(q)
    
    i <- i + 1
    
    q_list <- q_list[q]
    
  }
  #removing singles
  out <- out[which(as.vector(sapply(out, "length")) >= 2)]
  
  return(out)
   
}

put_na <- function(bool_set){
  if(bool_set){
    return(0)
  }
  else{
    return(NA)
  }
}

pull_queried <- function(gene_ids){
  
  whole_set <- list()
  for(i in 1:length(diffList)){
    #picking apart valid and invalid ids
    valid_ids <- which(gene_ids %in% diffList[[i]]$Symbol)
    
    nvalid_ids <- which(!(gene_ids %in% diffList[[i]]$Symbol))
    
    bool_set <- gene_ids %in% diffList[[i]]$Symbol
    #adding NAs to unavailable ids
    log2fc <- sapply(bool_set, FUN = put_na)
    
    #dealing with duplicate entries by taking the median from the database
    whole_ids <- diffList[[i]]$Symbol[diffList[[i]]$Symbol %in% gene_ids]
    dups <- make_dup_list(whole_ids)
    
    #browser()
    med_dups <- median_dups(logfc2 = diffList[[i]]$Log2FC, dupList = dups)
    
    #the duplicate location within the input id list
    dup_loc <- which(gene_ids %in% names(dups))
    
    log2fc[dup_loc] <- unlist(unname(med_dups))
    
    #removing valid_ids which are duplicates
    valid_ids <- valid_ids[which(!( valid_ids %in% dup_loc))]
    
    
    log2fc[valid_ids] <-  diffList[[i]]$Log2FC[which(diffList[[i]]$Symbol %in% gene_ids[valid_ids])]
    names(log2fc) <- gene_ids
    
    whole_set[names(diffList)[[i]]] <- list(log2fc)
    
    
    
  }
  #browser()
  dfs <- lapply(whole_set, data.frame, stringsAsFactors = FALSE)
  binded.ids <- bind_cols(dfs)
  rownames(binded.ids) <- gene_ids
  colnames(binded.ids) <- names(diffList)
  return(binded.ids)
  
}

logfcs <- pull_queried(query(varlist = list(pval = 0.001), smodel = "Random Effects Modeling"))
pca <- prcomp(na.omit(logfcs))



