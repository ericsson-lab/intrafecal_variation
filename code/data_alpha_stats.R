{
  source('code/load_data.R')
  
  library(vegan)
  library(microbiome)
}

alpha_stats <- table %>% 
  pivot_longer(-featureid, 
               names_to = "sampleid",
               values_to = "count") %>%
  group_by(sampleid) %>% 
  summarize(
    chao1 = as.double(richness(count, index = "chao1")),
    obs_richness = as.double(richness(count, index = "observed")),
    shannon = vegan::diversity(count, index = "shannon", base = exp(1)),
    simpson = vegan::diversity(count, index = "simpson")
  ) %>% 
  left_join(., metadata, by = 'sampleid')
