source('code/data_alpha_stats.R') # Load in alpha diversity stats
source('code/data_beta_diversity.R') # Load in beta diversity functions
library(rstatix)
library(lmerTest)
library(lme4)
library(ggprism)
library(usedist)

set.seed(1851)

theme_set(theme_prism())

# Filter to only GM Low/High samples
alpha_stats_gm <- alpha_stats %>% 
  filter(gm %in% c('Low', 'High')) %>% 
  mutate(gm = factor(gm, levels = c('Low', 'High')))

plot_alpha <- function(stat, ylab){
alpha_stats_gm %>% 
  ggplot(aes(x = gm, y = {{stat}}, shape = sex, color = gm)) +
  geom_point(size = 3,
             position = position_jitterdodge(jitter.width = 0.6,
                                             dodge.width = 0.9)) +
  stat_summary(geom = 'bar', fun = 'mean', fill = NA, color = 'black', linewidth = 1,
               width = 0.7, position = position_dodge(0.9)) +
  stat_summary(geom = 'errorbar', color = 'black', width = 0.2, linewidth = 1,
               fun.min = function(x) ifelse(mean(x) - sd(x) < 0, 0, mean(x) - sd(x)),
               fun.max = function(x) mean(x) + sd(x),
               position = position_dodge(0.9)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = c(expression(bold("GM"['Low'])), 
                              expression(bold('GM'['High'])))) +
  scale_color_manual(values = c('red', 'dodgerblue')) +
  scale_shape_manual(values = c(17, 16), name = 'Sex') +
  ylab(ylab) +
  theme_prism() +
  theme(
    legend.text = element_text(face = 'bold'),
    legend.title = element_text(face = 'bold', hjust = 0),
    axis.text = element_text(face = 'bold'),
    axis.title.y = element_text(face = 'bold'),
    axis.title.x = element_blank(),
    aspect.ratio = 3/1.75
  ) +
  guides(
    color = 'none',
    shape = guide_legend(override.aes = list(size = 4))
  )
}

# Plot alpha diversity figures
# Run Fit LMER and run ANOVA
plot_alpha(chao1, 'Chao1 Index')
ggsave('plots/chao1_overall.png', width = 5, height = 4)
chao1_lm <- alpha_stats_gm %>% 
  lmer(chao1 ~ gm * sex + (1|mouse/timepoint), data = .)
anova(chao1_lm)


plot_alpha(shannon, 'Shannon Index')
ggsave('plots/shannon_overall.png', width = 5, height = 4)
shannon_lm <- alpha_stats_gm %>% 
  lmer(shannon ~ gm * sex + (1|mouse/timepoint), data = .)
anova(shannon_lm)


plot_alpha(obs_richness, 'Observed Richness')
ggsave('plots/obs_richness_overall.png', width = 5, height = 4)
obs_richness_lm <- alpha_stats_gm %>% 
  lmer(obs_richness ~ gm * sex + (1|mouse/timepoint), data = .)
anova(obs_richness_lm)


plot_alpha(simpson, 'Simpson Index')
ggsave('plots/simpson_overall.png', width = 5, height = 4)
simpson_lm <- alpha_stats_gm %>% 
  lmer(simpson ~ gm * sex + (1|mouse/timepoint), data = .)
anova(simpson_lm)

# Generate transformed distance matrices
bc_dist <- generate_dist(table, 'bray')
j_dist <- generate_dist(table, 'jaccard')

gm_samples <- metadata %>% 
  filter(gm %in% c('Low', 'High'))

gm_bc_dist <-dist_subset(bc_dist, gm_samples$sampleid)
pcoa_bc <- generate_pcoa(gm_bc_dist)

gm_j_dist <-dist_subset(j_dist, gm_samples$sampleid)
pcoa_j <- generate_pcoa(gm_j_dist)

pcoa_bc[[1]] %>% 
  ggplot(aes(x = PCo1, y = PCo2, 
             color = factor(gm, levels = c('Low', 'High')),
             shape = sex)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('red', 'dodgerblue'),
                     labels = c(expression(bold('GM'['Low'])),
                                expression(bold('GM'['High']))),
                     name = 'GM') +
  scale_shape_manual(values = c(17, 16),
                     name = 'Sex') +
  xlab(glue::glue('PCo1 - {round(pcoa_bc[[2]][1],2)}%')) +
  ylab(glue::glue('PCo2 - {round(pcoa_bc[[2]][2],2)}%')) +
  theme(
    aspect.ratio = 1,
    legend.text = element_text(face = 'bold', size = 14),
    legend.title = element_text(face = 'bold', size = 14)
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 4)),
    color = guide_legend(override.aes = list(size = 4))
  )
ggsave('plots/bc_pcoa_overall.png', width = 6, height = 4)

# Create non-transformed distance matrix
bc_dist_adonis <-  table %>% 
  column_to_rownames(var = "featureid") %>% 
  t() %>% 
  vegdist(., method = 'bray')

gm_bc_dist_adonis <-dist_subset(bc_dist_adonis, gm_samples$sampleid)

# Set strata to only individual samples
perm = how(plots = with(pcoa_bc[[1]], Plots(strata = mouse:timepoint, type = 'free')),
           blocks = with(pcoa_bc[[1]], gm),
           nperm = 9999)

adonis2(gm_bc_dist_adonis ~ gm * sex, data = pcoa_bc[[1]], permutations = perm)


# Calculate coefficient of variance
sample <- alpha_stats_gm %>% 
  group_by(gm, mouse, timepoint) %>% 
  summarise(sample_chao1 = sd(chao1) / mean(chao1),
            sample_obs = sd(obs_richness) / mean(obs_richness),
            sample_shannon = sd(shannon) / mean(shannon),
            sample_simpson = sd(simpson) / mean(simpson),
  ) %>% 
  unite(mouse_timepoint, c('mouse', 'timepoint')) %>% 
  pivot_longer(-c(mouse_timepoint,gm )) %>% 
  select(-mouse_timepoint) %>% 
  separate(name, into = c('level', 'metric'))

mouse <- alpha_stats_gm %>% 
  group_by(gm, mouse) %>% 
  summarise(mouse_chao1 = sd(chao1) / mean(chao1),
            mouse_obs = sd(obs_richness) / mean(obs_richness),
            mouse_shannon = sd(shannon) / mean(shannon),
            mouse_simpson = sd(simpson) / mean(simpson),
  ) %>% 
  pivot_longer(-c(gm, mouse)) %>% 
  select(-mouse) %>% 
  separate(name, into = c('level', 'metric'))


cage <- alpha_stats_gm %>% 
  group_by(gm, cage) %>% 
  summarise(cage_chao1 = sd(chao1) / mean(chao1),
            cage_obs = sd(obs_richness) / mean(obs_richness),
            cage_shannon = sd(shannon) / mean(shannon),
            cage_simpson = sd(simpson) / mean(simpson),
  ) %>% 
  pivot_longer(-c(gm,cage)) %>% 
  select(-cage) %>% 
  separate(name, into = c('level', 'metric'))

# Bind CV data frames
cv_data <- sample %>% 
  rbind(., mouse) %>% 
  rbind(., cage)

cv_data %>% 
  mutate(level = str_to_sentence(level),
         level = case_match(level, 
                            'Sample' ~ 'Replicate',
                            .default = level),
         level = factor(level, levels = c('Replicate', 'Mouse', 'Cage')),
         metric = case_match(metric,
                             'chao1' ~ 'Chao1 Index',
                             'obs' ~ 'Observed Richness',
                             'shannon' ~ 'Shannon Index',
                             'simpson' ~ 'Simpson Index'
         )) %>% 
  ggplot(aes(x = level, y = value, color = gm)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.75), size = 3)+
  stat_summary(aes(group = gm), fun = 'mean', geom = 'bar', fill = NA, 
               color = 'black', width = 0.4, linewidth = 1, 
               position = position_dodge(0.75)) +
  stat_summary(aes(group = gm),
               geom = 'errorbar', color = 'black', width = 0.15, linewidth = 1,
               fun.min = function(x) ifelse(mean(x) - sd(x) < 0, 0, mean(x) - sd(x)),
               fun.max = function(x) mean(x) + sd(x),
               position = position_dodge(0.75)) +
  scale_color_manual(values = c('red', 'dodgerblue'),
                     labels = c(expression(bold('GM'['Low'])),
                                expression(bold('GM'['High']))),
                     name = 'GM') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  ylab('Coefficient of Variation') +
  facet_wrap(~metric, nrow = 1) +
  theme(
    legend.title = element_text(face = 'bold', size = 18),
    legend.text = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    strip.text = element_text(size = 16)
  )
ggsave('plots/cv_by_gm.png', width = 13, height = 4)

# run t test between GMs within each metric and level
cv_data %>% 
  group_by(level, metric) %>% 
  t_test(value ~ gm) %>% 
  arrange(p)

# Generate jaccard distance matrix for effect size calculations
j_dist_adonis <-  table %>% 
  column_to_rownames(var = "featureid") %>% 
  t() %>% 
  vegdist(., method = 'jaccard')

gm_j_dist_adonis <-dist_subset(j_dist_adonis, gm_samples$sampleid)

## Effect Size Calcualtions
bc_permanova <- adonis2(gm_bc_dist_adonis ~ gm / cage / mouse / timepoint, data = pcoa_bc[[1]],
        permutations = 9999)

## SOS-residual = 0.304
bc_permanova_es <- bc_permanova %>% 
  as_tibble(rownames = 'effect') %>% 
  drop_na %>% 
  mutate(es = SumOfSqs / (sum(SumOfSqs) + 0.304),
         value = 'PERMANOVA\n(Bray-Curtis)') %>% 
  rename(
    Effect = 'effect',
    p = `Pr(>F)`,
    es_manual = 'es'
  ) %>% 
  select(value, Effect, F, p, es_manual) 


j_permanova <- adonis2(gm_j_dist_adonis ~ gm / cage / mouse / timepoint, data = pcoa_bc[[1]],
                       permutations = 9999)
  ## SOS-residual = 1.087
j_permanova_es <- j_permanova %>% 
  as_tibble(rownames = 'effect') %>% 
  drop_na() %>% 
  mutate(es = SumOfSqs / (sum(SumOfSqs) + 1.087),
         value = 'PERMANOVA\n(Jaccard)') %>% 
  rename(
    Effect = 'effect',
    p = `Pr(>F)`,
    es_manual = 'es'
  ) %>% 
  select(value, Effect, F, p, es_manual)


### Performing nested ANOVA and manually calculating eta squared due to
### error in rstatix code: https://github.com/kassambara/rstatix/issues/132
chao1 <- alpha_stats_gm %>% 
  anova_test(chao1 ~  gm / cage / mouse / timepoint, effect.size = 'pes',
             detailed = T) %>% 
  add_significance() %>% 
  as_tibble() %>% 
  mutate(value = 'Chao1') %>% 
  mutate(pes_manual = SSn / (SSn + SSd),
         es_manual = SSn / (sum(SSn) + SSd)) %>% 
  select(value, Effect, F, p, es_manual)


obs_richness <- alpha_stats_gm %>% 
  anova_test(obs_richness ~ gm / cage/ mouse / timepoint, effect.size = 'pes',
             detailed = T) %>% 
  add_significance() %>% 
  as_tibble() %>% 
  mutate(value = 'Observed\nRichness') %>% 
  mutate(pes_manual = SSn / (SSn + SSd),
         es_manual = SSn / (sum(SSn) + SSd)) %>% 
  select(value, Effect, F, p, es_manual)


shannon <- alpha_stats_gm %>% 
  anova_test(shannon ~ gm / cage / mouse / timepoint, effect.size = 'pes',
             detailed = T, ) %>% 
  add_significance() %>% 
  as_tibble() %>% 
  mutate(value = 'Shannon') %>% 
  mutate(pes_manual = SSn / (SSn + SSd),
         es_manual = SSn / (sum(SSn) + SSd)) %>% 
  select(value, Effect, F, p, es_manual)


simpson <- alpha_stats_gm %>% 
  anova_test(simpson ~ gm /cage / mouse / timepoint, effect.size = 'pes',
             detailed = T, type = 1) %>% 
  add_significance() %>% 
  as_tibble() %>% 
  mutate(value = 'Simpson') %>% 
  mutate(pes_manual = SSn / (SSn + SSd),
         es_manual = SSn / (sum(SSn) + SSd)) %>% 
  select(value, Effect, F, p, es_manual)


effectSize <- rbind(chao1, shannon, simpson, obs_richness,
                    bc_permanova_es, j_permanova_es) %>% 
  mutate(Effect = case_match(Effect,
                             'gm' ~ 'GM',
                             'gm:cage' ~ 'Cage',
                             'gm:cage:mouse' ~ 'Mouse',
                             'gm:cage:mouse:timepoint' ~ 'Replicate'),
         Effect = factor(Effect, levels = c('GM', 'Cage', 'Mouse', 'Replicate')),
         value = factor(value, levels = c('Observed\nRichness', 'Chao1', 'Shannon', 'Simpson',
                                          'PERMANOVA\n(Bray-Curtis)', 'PERMANOVA\n(Jaccard)'))) %>% 
  add_significance(p.col = 'p', output.col = 'p.signif',
                   cutpoints = c(0,0.001, 0.01, 0.05, 1),
                   symbols = c('***', '**', '*', NA))

effectSize %>% 
  ggplot(aes(x = value, y = fct_rev(Effect), fill = es_manual)) +
  geom_tile(color = 'white', width=1, height = 1) +
  geom_text(aes(label = p.signif), size = 10, color = 'white',
            position = position_nudge(y = -0.1)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_gradient(low = '#c4a8e5', high = '#36135f', limits = c(0, 1),
                      breaks = c(0, 0.5, 1), name = 'Eta\nSquared') +
  theme(
    axis.title = element_blank(),
    legend.text = element_text(face = 'bold'),
    legend.title = element_text(face = 'bold', hjust = 0),
    aspect.ratio = 4/6 ,
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1)
  ) +
  guides(
    fill = guide_colorbar(ticks.colour = NA)
  )
ggsave('plots/pes_both_gms.png', width = 8, height = 5)

effectSize %>% 
  group_by(Effect) %>% 
  summarize(es_mean = mean(es_manual),
            es_sd = sd(es_manual)) %>% 
  clipr::write_clip()



# Effect size calculations by GM
chao1_gm <- alpha_stats_gm %>%
  group_by(gm) %>% 
  anova_test(chao1 ~  cage / mouse / timepoint, effect.size = 'pes',
             detailed = T) %>% 
  add_significance() %>% 
  as_tibble() %>% 
  mutate(value = 'Chao1') %>% 
  mutate(pes_manual = SSn / (SSn + SSd),
         es_manual = SSn / (sum(SSn) + SSd)) %>% 
  select(gm, value, Effect, F, p, es_manual)

obs_richness_gm <- alpha_stats_gm %>% 
  group_by(gm) %>% 
  anova_test(obs_richness ~ cage/ mouse / timepoint, effect.size = 'pes',
             detailed = T) %>% 
  add_significance() %>% 
  as_tibble() %>% 
  mutate(value = 'Observed\nRichness') %>% 
  mutate(pes_manual = SSn / (SSn + SSd),
         es_manual = SSn / (sum(SSn) + SSd)) %>% 
  select(gm, value, Effect, F, p, es_manual)


shannon_gm <- alpha_stats_gm %>%
  group_by(gm) %>% 
  anova_test(shannon ~ cage / mouse / timepoint, effect.size = 'pes',
             detailed = T, ) %>% 
  add_significance() %>% 
  as_tibble() %>% 
  mutate(value = 'Shannon') %>% 
  mutate(pes_manual = SSn / (SSn + SSd),
         es_manual = SSn / (sum(SSn) + SSd)) %>% 
  select(gm, value, Effect, F, p, es_manual)


simpson_gm <- alpha_stats_gm %>% 
  group_by(gm) %>% 
  anova_test(simpson ~ cage / mouse / timepoint, effect.size = 'pes',
             detailed = T, type = 1) %>% 
  add_significance() %>% 
  as_tibble() %>% 
  mutate(value = 'Simpson') %>% 
  mutate(pes_manual = SSn / (SSn + SSd),
         es_manual = SSn / (sum(SSn) + SSd)) %>% 
  select(gm, value, Effect, F, p, es_manual)



gm_low_samples_list <- gm_samples %>% 
  filter(gm == 'Low')
gm_high_samples_list <- gm_samples %>% 
  filter(gm == 'High')


## GM Low
#  Jaccard
gm_low_j_dist <-dist_subset(j_dist, gm_low_samples_list$sampleid)
pcoa_low_j <- generate_pcoa(gm_low_j_dist)

gm_low_j_permanova <- adonis2(gm_low_j_dist ~ cage / mouse / timepoint, data = pcoa_low_j[[1]],
                        permutations = 9999)
## SOS-residual = 0.5157
gm_low_j_permanova_es <- gm_low_j_permanova %>% 
  as_tibble(rownames = 'effect') %>% 
  drop_na %>% 
  mutate(es = SumOfSqs / (sum(SumOfSqs) + 0.5157),
         value = 'PERMANOVA\n(Jaccard)',
         gm = "Low") %>% 
  rename(
    Effect = 'effect',
    p = `Pr(>F)`,
    es_manual = 'es'
  ) %>% 
  select(gm, value,  Effect, F, p, es_manual) 

#BC
gm_low_bc_dist <-dist_subset(bc_dist_adonis, gm_low_samples_list$sampleid)
pcoa_low_bc <- generate_pcoa(gm_low_bc_dist)

gm_low_bc_permanova <- adonis2(gm_low_bc_dist ~ cage / mouse / timepoint, data = pcoa_low_bc[[1]],
                              permutations = 9999)
## SOS-residual = 0.1439
gm_low_bc_permanova_es <- gm_low_bc_permanova %>% 
  as_tibble(rownames = 'effect') %>% 
  drop_na %>% 
  mutate(es = SumOfSqs / (sum(SumOfSqs) + 0.1439),
         value = 'PERMANOVA\n(Bray-Curtis)',
         gm = 'Low') %>% 
  rename(
    Effect = 'effect',
    p = `Pr(>F)`,
    es_manual = 'es'
  ) %>% 
  select(gm, value, Effect, F, p, es_manual) 


## GM High
gm_high_j_dist <-dist_subset(j_dist_adonis, gm_high_samples_list$sampleid)
pcoa_high_j <- generate_pcoa(gm_high_j_dist)

gm_high_j_permanova <- adonis2(gm_high_j_dist ~ cage / mouse / timepoint, data = pcoa_high_j[[1]],
                              permutations = 9999)
## SOS-residual = 0.5718
gm_high_j_permanova_es <- gm_high_j_permanova %>% 
  as_tibble(rownames = 'effect') %>% 
  drop_na %>% 
  mutate(es = SumOfSqs / (sum(SumOfSqs) + 0.5718),
         value = 'PERMANOVA\n(Jaccard)',
         gm = 'High') %>% 
  rename(
    Effect = 'effect',
    p = `Pr(>F)`,
    es_manual = 'es'
  ) %>% 
  select(gm, value, Effect, F, p, es_manual) 

#BC
gm_high_bc_dist <-dist_subset(bc_dist_adonis, gm_high_samples_list$sampleid)
pcoa_high_bc <- generate_pcoa(gm_high_bc_dist)

gm_high_bc_permanova <- adonis2(gm_high_bc_dist ~ cage / mouse / timepoint, data = pcoa_high_bc[[1]],
                               permutations = 9999)
## SOS-residual = 0.1604
gm_high_bc_permanova_es <- gm_high_bc_permanova %>% 
  as_tibble(rownames = 'effect') %>% 
  drop_na %>% 
  mutate(es = SumOfSqs / (sum(SumOfSqs) + 0.1604),
         value = 'PERMANOVA\n(Bray-Curtis)',
         gm = 'High') %>% 
  rename(
    Effect = 'effect',
    p = `Pr(>F)`,
    es_manual = 'es'
  ) %>% 
  select(gm, value, Effect, F, p, es_manual)

effectSize_gm <- rbind(chao1_gm, shannon_gm, simpson_gm, obs_richness_gm,
                       gm_low_j_permanova_es, gm_low_bc_permanova_es,
                       gm_high_j_permanova_es, gm_high_bc_permanova_es) %>% 
  mutate(Effect = case_match(Effect,
                             'cage' ~ 'Cage',
                             'cage:mouse' ~ 'Mouse',
                             'cage:mouse:timepoint' ~ 'Replicate'),
         Effect = factor(Effect, levels = c('Cage', 'Mouse', 'Replicate')),
         value = factor(value, levels = c('Observed\nRichness', 'Chao1', 'Shannon', 'Simpson',
                                          'PERMANOVA\n(Bray-Curtis)', 'PERMANOVA\n(Jaccard)'))) %>% 
  add_significance(p.col = 'p', output.col = 'p.signif',
                   cutpoints = c(0,0.001, 0.01, 0.05, 1),
                   symbols = c('***', '**', '*', NA)) 
effectSize_gm %>% 
  arrange(-es_manual)

plot_gm_eta_sq <- function(gm, low_col, high_col){
  effectSize_gm %>% 
    filter(gm == {{gm}}) %>% 
    ggplot(aes(x = value, y = fct_rev(Effect), fill = es_manual),) +
    geom_tile(color = 'white', width=1, height = 1) +
    geom_text(aes(label = p.signif), size = 10, color = 'white',
              position = position_nudge(y = -0.1)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0)) +
    scale_fill_gradient(low = low_col, high = high_col, limits = c(0, 0.75),
                        breaks = c(0, 0.25, 0.5, 0.75), name = 'Eta\nSquared') +
    theme(
      axis.title = element_blank(),
      legend.text = element_text(face = 'bold'),
      legend.title = element_text(face = 'bold', hjust = 0),
      aspect.ratio = 3/6,
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1)
    ) +
    guides(
      fill = guide_colorbar(ticks.colour = NA)
    )
}
  
plot_gm_eta_sq(gm = 'High', low_col = '#c1d0f6', high_col = '#0f2a72')
ggsave('plots/es_gm_high.png', width = 8, height = 5)
plot_gm_eta_sq(gm = 'Low', low_col = '#FFa5a1', high_col = '#9b0700')
ggsave('plots/es_gm_low.png', width = 8, height = 5)

effectSize_gm %>% 
  group_by(gm, Effect) %>% 
  summarize(
    mean_es = mean(es_manual),
    sd_es = sd(es_manual)
  ) %>% 
  clipr::write_clip()


## Measuring precision of beta diversity
sample_metadata <- gm_samples %>% 
  filter(gm %in% c('Low', 'High')) %>% 
  unite(col = 'mouse_timepoint', c('mouse', 'timepoint'), remove = F)

sample_list <- sample_metadata %>% 
  pull(mouse_timepoint) %>% 
  unique()
length(sample_list)

result = tibble()

for (i in 1:length(sample_list)) {
  technical_samples <- sample_metadata %>% 
    filter(mouse_timepoint == sample_list[i]) 
  
  technical_dist <- dist_subset(bc_dist, idx = technical_samples$sampleid)
  result_add = tibble(mouse_timepoint = sample_list[i], avg_dist = mean(technical_dist))
  result = rbind(result, result_add)
}

replicate_result <- result %>% 
  left_join(sample_metadata %>% 
              select(mouse_timepoint, gm, sex) %>% 
              distinct()) %>% 
  select(gm, avg_dist, mouse_timepoint) %>% 
  rename(identifier = 'mouse_timepoint') %>% 
  mutate(level = 'replicate')


mouse_list <- sample_metadata %>% 
  pull(mouse) %>% 
  unique()
length(mouse_list)

result = tibble()

for (i in 1:length(mouse_list)) {
  technical_samples <- sample_metadata %>% 
    filter(mouse == mouse_list[i]) 
  
  technical_dist <- dist_subset(bc_dist, idx = technical_samples$sampleid)
  result_add = tibble(mouse = mouse_list[i], avg_dist = mean(technical_dist))
  result = rbind(result, result_add)
}

mouse_result <- result %>% 
  left_join(sample_metadata %>% 
              select(mouse, gm, sex) %>% 
              distinct()) %>% 
  select(gm, avg_dist, mouse) %>% 
  rename(identifier = 'mouse') %>% 
  mutate(level = 'mouse')


cage_list <- sample_metadata %>% 
  pull(cage) %>% 
  unique()
length(cage_list)

result = tibble()

for (i in 1:length(cage_list)) {
  technical_samples <- sample_metadata %>% 
    filter(cage == cage_list[i]) 
  
  technical_dist <- dist_subset(bc_dist, idx = technical_samples$sampleid)
  result_add = tibble(cage = cage_list[i], avg_dist = mean(technical_dist))
  result = rbind(result, result_add)
}

cage_result <- result %>% 
  left_join(sample_metadata %>% 
              select(cage, gm, sex) %>% 
              distinct()) %>% 
  select(gm, avg_dist, cage) %>% 
  rename(identifier = 'cage') %>% 
  mutate(level = 'cage')


gmLow_samples <- sample_metadata %>% 
  filter(gm == "Low")

gmLow_dist <- dist_subset(bc_dist, idx = gmLow_samples$sampleid)
gmLow_result = tibble(gm = 'Low', avg_dist = mean(gmLow_dist), sd = sd(gmLow_dist),
                      identifier = 'Low', level = 'Low')

gmHigh_samples <- sample_metadata %>% 
  filter(gm == "High")
gmHigh_dist <- dist_subset(bc_dist, idx = gmHigh_samples$sampleid)
gmHigh_result = tibble(gm = 'High', avg_dist = mean(gmHigh_dist),sd = sd(gmHigh_dist),
                       identifier = 'High', level = 'High')


# Combine beta diversity measurements
bc_intralevel_distances <- replicate_result %>% 
  rbind(mouse_result) %>% 
  rbind(cage_result) %>% 
  mutate(level = factor(level, 
                        levels = c('replicate', 'mouse', 'cage', 'Low', 'High')),
         gm = factor(gm, levels = c('Low', 'High'))) %>% 
  rename(bc_avg_dist = 'avg_dist')

result = tibble()

for (i in 1:length(sample_list)) {
  technical_samples <- sample_metadata %>% 
    filter(mouse_timepoint == sample_list[i]) 
  
  technical_dist <- dist_subset(j_dist, idx = technical_samples$sampleid)
  result_add = tibble(mouse_timepoint = sample_list[i], avg_dist = mean(technical_dist))
  result = rbind(result, result_add)
}

replicate_result <- result %>% 
  left_join(sample_metadata %>% 
              select(mouse_timepoint, gm, sex) %>% 
              distinct()) %>% 
  select(gm, avg_dist, mouse_timepoint) %>% 
  rename(identifier = 'mouse_timepoint') %>% 
  mutate(level = 'replicate')

result = tibble()

for (i in 1:length(mouse_list)) {
  technical_samples <- sample_metadata %>% 
    filter(mouse == mouse_list[i]) 
  
  technical_dist <- dist_subset(j_dist, idx = technical_samples$sampleid)
  result_add = tibble(mouse = mouse_list[i], avg_dist = mean(technical_dist))
  result = rbind(result, result_add)
}

mouse_result <- result %>% 
  left_join(sample_metadata %>% 
              select(mouse, gm, sex) %>% 
              distinct()) %>% 
  select(gm, avg_dist, mouse) %>% 
  rename(identifier = 'mouse') %>% 
  mutate(level = 'mouse')


result = tibble()

for (i in 1:length(cage_list)) {
  technical_samples <- sample_metadata %>% 
    filter(cage == cage_list[i]) 
  
  technical_dist <- dist_subset(j_dist, idx = technical_samples$sampleid)
  result_add = tibble(cage = cage_list[i], avg_dist = mean(technical_dist))
  result = rbind(result, result_add)
}

cage_result <- result %>% 
  left_join(sample_metadata %>% 
              select(cage, gm, sex) %>% 
              distinct()) %>% 
  select(gm, avg_dist, cage) %>% 
  rename(identifier = 'cage') %>% 
  mutate(level = 'cage')

j_intralevel_distances <- replicate_result %>% 
  rbind(mouse_result) %>% 
  rbind(cage_result) %>% 
  mutate(level = factor(level, 
                        levels = c('replicate', 'mouse', 'cage', 'Low', 'High')),
         gm = factor(gm, levels = c('Low', 'High'))) %>% 
  rename(j_avg_dist = 'avg_dist')

bc_intralevel_distances %>% 
  left_join(., j_intralevel_distances) %>% 
  mutate(level = str_to_sentence(level)) %>% 
  filter(!level %in% c('Low', 'High')) %>% 
  pivot_longer(-c(gm, identifier, level), values_to = 'avg_dist',
               names_to = 'dist') %>% 
  ggplot(aes(x = factor(level, levels = c('Replicate', 'Mouse', 'Cage')), 
             y = avg_dist, 
             color = gm))+
  geom_point(position = position_jitter( 0.2), 
             size = 3, color = '#966fd6')  +
  stat_summary(position = position_dodge( 0.75),
               geom = 'bar', fun = 'mean',
               fill = NA, color ='black',
               width = 0.6, linewidth = 1) +
  stat_summary(geom = 'errorbar',
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               color ='black', linewidth = 1,
               aes(group = dist), width = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0,0.01)),
                     limits = c(0, 0.5)) +
  ylab('Dissimilarity') +
  theme(
    aspect.ratio = 1,
    axis.title.x = element_blank(),
    legend.position = 'none'
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) +
  facet_wrap(~dist)

ggsave('plots/j_bc_distances.png', width = 6, height = 4)

bc_intralevel_distances %>% 
  left_join(., j_intralevel_distances) %>% 
  group_by(level) %>% 
  summarise(mean_bc = mean(bc_avg_dist),
            sd_bc = sd(bc_avg_dist),
            mean_j = mean(j_avg_dist),
            sd_j = sd(j_avg_dist)) %>% 
  clipr::write_clip()


bc_intralevel_distances %>% 
  left_join(., j_intralevel_distances) %>% 
  mutate(level = str_to_sentence(level)) %>% 
  filter(!level %in% c('Low', 'High')) %>% 
  pivot_longer(-c(gm, identifier, level), values_to = 'avg_dist',
               names_to = 'dist') %>% 
  ggplot(aes(x = factor(level, levels = c('Replicate', 'Mouse', 'Cage')), 
             y = avg_dist, 
             color = gm)) +
  geom_point(position = position_jitter(0.2), 
             size = 3) + 
  stat_summary(aes(group = dist),
               geom = 'bar', fun = 'mean',
               fill = NA, color ='black',
               width = 0.6, linewidth = 1)+
  stat_summary(geom = 'errorbar',
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               color ='black', linewidth = 1,
               aes(group = dist), width = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0,0.01)),
                     limits = c(0, 0.5)) +
  scale_shape_manual(values = c(15, 17), 
                     name = 'Distance Metric',
                     labels = c('Bray-Curtis',
                                'Jaccard')) +
  ylab('Dissimilarity') +
  scale_color_manual(values = c('red', 'dodgerblue'),
                     labels = c(expression(bold('GM'['Low'])),
                                expression(bold('GM'['High']))),
                     name = 'GM') +
  theme(
    aspect.ratio = 1,
    axis.title.x = element_blank(),
    legend.position = 'none',
    strip.text = element_blank()
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) +
  facet_wrap(dist~gm, scales = 'free_x')

ggsave('plots/bc_distances_by_gm.png', width = 10, height = 6)


bc_intralevel_distances %>% 
  left_join(., j_intralevel_distances) %>% 
  group_by(level, gm) %>% 
  summarize(mean_bc = mean(bc_avg_dist),
            sd_bc = sd(bc_avg_dist),
            mean_j = mean(j_avg_dist),
            sd_j = sd(j_avg_dist),) %>% 
  drop_na() %>% 
  clipr::write_clip()



# Sample Size Calculations
# Need to clean up the output function
plot_effectsize <- function(stat, diff_set){
  
  sample_list = tibble(n_samples = rep(c(1:10),3),
                       diff = c(rep(diff_set[1], 10), rep(diff_set[2], 10), rep(diff_set[3], 10)))
 
   # replicate (technical)
  technical_var <- alpha_stats_gm %>% 
    select(sampleid, gm, cage, mouse, timepoint, {{stat}}) %>% 
    group_by(mouse, timepoint) %>% 
    mutate(mean = mean({{stat}}),
           mean_dif = {{stat}} - mean,
           squared_mean_diff = mean_dif^2) %>% 
    ungroup() %>% 
    summarize(var = sum(squared_mean_diff) / (length(squared_mean_diff) - 1)) %>% 
    pull(var)
  # Sample
  sample_var <- alpha_stats_gm %>% 
    select(sampleid, gm, cage, mouse, timepoint, {{stat}}) %>% 
    group_by(mouse, timepoint) %>% 
    summarise(stat = mean({{stat}}),.groups = 'drop') %>% 
    group_by(mouse) %>%
    mutate(mean = mean(stat),
           mean_dif = stat - mean,
           squared_mean_diff = mean_dif^2) %>% 
    ungroup() %>% 
    summarize(var = sum(squared_mean_diff) / (length(squared_mean_diff) - 1)) %>% 
    pull(var)
  
  # Mouse
  mouse_var <- alpha_stats_gm %>% 
    select(sampleid, gm, cage, mouse, timepoint, {{stat}}) %>% 
    group_by(mouse) %>% 
    summarise(stat = mean({{stat}}),.groups = 'drop') %>% 
    mutate(mean = mean(stat),
           mean_dif = stat - mean,
           squared_mean_diff = mean_dif^2)  %>% 
    ungroup() %>% 
    summarize(var = sum(squared_mean_diff) / (length(squared_mean_diff) - 1)) %>% 
    pull(var)
  
  n_tech_replicate  = 1
  
  power_list <- sample_list %>% 
    mutate(var = mouse_var + (sample_var / n_samples) + 
             (technical_var / (n_samples * n_tech_replicate)),
           d = diff/sqrt(var)
           )
  power_list %>% 
    ggplot(aes(x = (n_samples), y = d, group = diff, linetype = as.factor(diff))) +
    geom_line(linewidth = 1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(breaks = c(seq(1,10,1))) +
    scale_linetype(name = 'Difference') +
    xlab('Samples Per Mouse') +
    ylab("Cohen's (*d*)") +
    theme(
      aspect.ratio = 2/1,
      axis.title.y = ggtext::element_markdown()
    )
  return(power_list)
}

tmp1 <- plot_effectsize(stat = chao1,diff_set = c(96, 50, 25)) 

sample_size <- array()

for (i in 1:length(tmp1$d)) {
  tmp <- pwr::pwr.t.test(d = tmp1$d[i], sig.level = 0.05, power = 0.8, alternative = 'two.sided')$n 
  sample_size[[i]] <- tmp 
}

# Perform cost analysis
tmp2 <- tmp1 %>% 
  cbind(sample_size) %>% 
  mutate(cost = 50 * n_samples * sample_size) 
  
tmp2 %>% 
  ggplot(aes(x = n_samples, y = cost, linetype = as.factor(diff))) +
  geom_line(linewidth = 2) +
  scale_y_continuous(labels = scales::dollar, 
                     limits = c(0, 4e4),
                     expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(breaks = 1:10) +
  scale_linetype_manual(values = c(3,2,1)) +
  xlab('Samples / Mouse') +
  ylab("Cost / Group") +
  theme(
    aspect.ratio = 1,
    legend.position = 'none'
  )
ggsave('plots/cost.png', width = 4.25, height = 4.25)

tmp2 %>% 
  ggplot(aes(x = n_samples, y = sample_size, linetype = as.factor(diff))) +
  geom_line(linewidth = 2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     limits = c(0, 80)) +
  scale_x_continuous(breaks = 1:10) +
  scale_linetype_manual(values = c(3,2,1)) +
  xlab('Samples / Mouse') +
  ylab("Mice / Group") +
  theme(
    aspect.ratio = 1,
    legend.position = 'none'
  )
ggsave('plots/sample_size.png', width = 4, height = 4)

tmp2 %>% 
  ggplot(aes(x = n_samples, y = d, linetype = as.factor(diff))) +
  geom_line(linewidth = 2.1) +
  scale_y_continuous(limits = c(0, 2),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_linetype_manual(values = c(3,2,1)) +
  scale_x_continuous(breaks = 1:10) +
  xlab('Samples / Mouse') +
  ylab("Cohen's d") + 
  theme(
    aspect.ratio = 1,
    legend.position = 'none'
  )
ggsave('plots/cohens_d.png', width = 4, height = 4)

tmp2 %>% 
  ggplot(aes(x = n_samples, y = cost, linetype = factor(diff))) +
  geom_line(linewidth = 1) +
  scale_linetype_manual(values = c(3,2,1), name = 'Difference in\nMeans') +
  guides(
    linetype = guide_legend(override.aes = list(linewidth = 0.75))
  ) +
  theme(
    legend.title = element_text(face = 'bold', size = 10),
    legend.text = element_text(face = 'bold', size = 10)
  )
ggsave('plots/legend.png', width = 8, height = 8)


tmp2 %>% 
  group_by(diff) %>% 
  mutate(
    d_percent_increase = (d/first(d)-1),
    d_sample_size_decrease = (sample_size/first(sample_size)-1),
    d_cost_increase = (cost/first(cost)-1)
  ) %>% 
  write_tsv('stats/effectSize_sampleSize_cost.tsv')


# Sample size calcualtions for alpha diveristy stats
plot_sample_size <- function(stat, diff_range, metric){

  sample_list = tibble(diff = diff_range)
  
  # replicate (technical)
  technical_var <- alpha_stats_gm %>% 
    select(sampleid, gm, cage, mouse, timepoint, {{stat}}) %>% 
    group_by(mouse, timepoint) %>% 
    mutate(mean = mean({{stat}}),
           mean_dif = {{stat}} - mean,
           squared_mean_diff = mean_dif^2) %>% 
    ungroup() %>% 
    summarize(var = sum(squared_mean_diff) / (length(squared_mean_diff) - 1)) %>% 
    pull(var)
  
  # Sample
  sample_var <- alpha_stats_gm %>% 
    select(sampleid, gm, cage, mouse, timepoint, {{stat}}) %>% 
    group_by(mouse, timepoint) %>% 
    summarise(stat = mean({{stat}}),.groups = 'drop') %>% 
    group_by(mouse) %>%
    mutate(mean = mean(stat),
           mean_dif = stat - mean,
           squared_mean_diff = mean_dif^2) %>% 
    ungroup() %>% 
    summarize(var = sum(squared_mean_diff) / (length(squared_mean_diff) - 1)) %>% 
    pull(var)
  
  # Mouse
  mouse_var <- alpha_stats_gm %>% 
    select(sampleid, gm, cage, mouse, timepoint, {{stat}}) %>% 
    group_by(mouse) %>% 
    summarise(stat = mean({{stat}}),.groups = 'drop') %>% 
    mutate(mean = mean(stat),
           mean_dif = stat - mean,
           squared_mean_diff = mean_dif^2)  %>% 
    ungroup() %>% 
    summarize(var = sum(squared_mean_diff) / (length(squared_mean_diff) - 1)) %>% 
    pull(var)
  
  n_samples = 1
  n_tech_replicate  = 1
  
  power_list <- sample_list %>% 
    mutate(var = mouse_var + (sample_var / n_samples) + 
             (technical_var / (n_samples * n_tech_replicate)),
           d = diff/sqrt(var)
    )
  
  sample_size <- array()

  for (i in 1:length(power_list$d)) {
    
    tmp <- pwr::pwr.t.test(d = power_list$d[i], 
                           sig.level = 0.05, 
                           power = 0.8, 
                           alternative = 'two.sided')$n 
    sample_size[[i]] <- tmp 
  }
  power_list <- power_list %>% 
    cbind(., sample_size)
  
  power_list %>% 
    ggplot(aes(x = diff, y = log10(sample_size))) +
    geom_line(linewidth = 1) +
    scale_linetype(name = 'Difference') +
    xlab(glue::glue('{metric}\nDifference in Means')) +
    ylab("Mice / Group<br>(Log<sub>10</sub><i>n</i>)") +
    theme(
      aspect.ratio = 1,
      axis.title.y = ggtext::element_markdown() 
    )

}

a_chao1 <- plot_sample_size(chao1, 1:250, 'Chao1 Index') + 
  scale_y_continuous(limits = c(0, 5),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(limits = c(0, 250),
                     expand = expansion(mult = c(0, 0.05))) 
  

b_obs_richness <- plot_sample_size(obs_richness, 1:250, 'Observed Richness') + 
  scale_y_continuous(limits = c(0, 5),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(limits = c(0, 250),
                     expand = expansion(mult = c(0, 0.05))) 

c_shannon <- plot_sample_size(shannon, seq(0.01, 1., 0.01), 'Shannon Index') +
  scale_y_continuous(limits = c(0, 5),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.05))) 

d_simpson <- plot_sample_size(simpson, seq(0.001, 0.05, 0.001), 'Simpson Index') +
  scale_y_continuous(limits = c(0, 5),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(limits = c(0, 0.05),
                     expand = expansion(mult = c(0, 0.05))) 




cowplot::plot_grid(a_chao1, b_obs_richness, c_shannon, d_simpson,
                   nrow = 2, ncol = 2) 
ggsave('plots/sample_size.png', width = 7, height = 7, bg = 'white')






