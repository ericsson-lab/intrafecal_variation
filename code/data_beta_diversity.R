{
  library(vegan)
  library(ape)
  source("code/load_data.R")
  
}

generate_dist <- function(table, distance = "bray"){
  
  table_t <- table %>% 
    column_to_rownames(var = "featureid") %>% 
    t()
  
  table.transformed <- table_t^(1/4)
  
  dist <- vegan::vegdist(table.transformed, method= distance)
}

generate_pcoa <- function(dist) {
  pcoa_res <- ape::pcoa(dist, correction = "cailliez") 
  
  p_var <- (pcoa_res$values$Eigenvalues/pcoa_res$trace)*100
  
  pcoa_vectors <- pcoa_res$vectors %>% 
    as_tibble(rownames = "sampleid") %>% 
    select(sampleid, Axis.1, Axis.2, Axis.3)
  
  colnames(pcoa_vectors) <- c("sampleid", "PCo1", "PCo2", "PCo3")
  
  pcoa.metadata <- left_join(pcoa_vectors, metadata, by = "sampleid") 

  output <- list(pcoa.metadata, p_var)
  return(output)
}





