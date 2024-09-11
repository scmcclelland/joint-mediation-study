# Functions for Analysis Script 
# file name:   analysis-functions.R
# Author:      Shelby C McClelland
# Created:     20 December 2020
# Last Update: 18 April 2024
# Description: This file contains functions to analyze data for McClelland and
#              Schipanski (XXXX). YYYY. Ecology Letters.
#-------------------------------------------------------------------------------
boxplot_check = function(dt, variable, level) {
  require(lattice)
  bw = bwplot(get(variable) ~ Trt|Depth + Year, data = dt)
  return(bw)
}
boxplot_check2 = function(dt, variable) {
  require(lattice)
  bw = bwplot(get(variable) ~ Trt|Year, data = dt)
  return(bw)
}
repeated_MM = function(dt, variable, level) {
  require(lme4)
  require(emmeans)
  print('Running lme model')
  model     = lmer(get(variable) ~ Trt*Year + (1|Block) 
                + (1|Sample_ID), data = dt[Depth %in% level])
  infl      = influence(model, obs  = TRUE)
  infl.cook = cooks.distance(model, infl)
  plot_m    = plot(model, which = "cook")
  print('Generating ANOVA')
  anova_r   = anova(model, ddf = 'Kenward-Roger')
  print('Generating emmeans')
  emmeans_r = emmeans(model, pairwise ~ Trt|Year)
  result    = list('model_description' = paste(variable, level, sep = '-'),
                'model_summary'    = summary(model), 
                'influence'        = infl, 
                'influence_cook'   = infl.cook,
                'model_plot'       = plot_m,
                'anova'            = anova_r,
                'emmeans'          = emmeans_r) 
  return(result)
}
repeated_MM2 = function(dt, variable) {
  require(lme4)
  require(emmeans)
  print('Running lme model')
  model     = lmer(get(variable) ~ Trt*Year + (1|Block) 
                   + (1|Sample_ID), data = dt)
  infl      = influence(model, obs  = TRUE)
  infl.cook = cooks.distance(model, infl)
  plot_m    = plot(model, which = "cook")
  print('Generating ANOVA')
  anova_r   = anova(model, ddf = 'Kenward-Roger')
  print('Generating emmeans')
  emmeans_r = emmeans(model, pairwise ~ Trt|Year)
  result    = list('model_description' = paste(variable, sep = '-'),
                   'model_summary'    = summary(model), 
                   'influence'        = infl, 
                   'influence_cook'   = infl.cook,
                   'model_plot'       = plot_m,
                   'anova'            = anova_r,
                   'emmeans'          = emmeans_r) 
  return(result)
}
repeated_MM3 = function(dt, variable) {
  require(lme4)
  require(emmeans)
  print('Running lme model')
  model     = lmer(get(variable) ~ Trt*Year*Substrate + (1|Block) 
                   + (1|Sample_ID), data = dt)
  infl      = influence(model, obs  = TRUE)
  infl.cook = cooks.distance(model, infl)
  plot_m    = plot(model, which = "cook")
  print('Generating ANOVA')
  anova_r   = anova(model, ddf = 'Kenward-Roger')
  print('Generating emmeans')
  emmeans_r = emmeans(model, pairwise ~ Trt|Substrate + Year)
  result    = list('model_description' = paste(variable, sep = '-'),
                   'model_summary'    = summary(model), 
                   'influence'        = infl, 
                   'influence_cook'   = infl.cook,
                   'model_plot'       = plot_m,
                   'anova'            = anova_r,
                   'emmeans'          = emmeans_r) 
  return(result)
}
MM = function(dt, variable, level) {
  require(lme4)
  require(emmeans)
  print('Running lme model')
  model     = lmer(get(variable) ~ Trt + (1|Block), data = dt[Depth %in% level])
  plot_m    = plot(model)
  print('Generating ANOVA')
  anova_r   = anova(model, ddf = 'Kenward-Roger')
  result    = list('model_description' = paste(variable, level, sep = '-'),
                   'model_summary' = summary(model), 
                   'model_plot'       = plot_m,
                   'anova'            = anova_r) 
  return(result)
}
