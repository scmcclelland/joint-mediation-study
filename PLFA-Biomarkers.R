# PLFA Biomarkers Analysis Script 
# Author: Shelby C McClelland
# Created:     20 December 2020
# Last Update: 22 February 2025
# Description: This file analyzes soil chemical property data.
#-------------------------------------------------------------------------------
## Analysis notes:
# Repeated measures analysis
# Block and plot treated as random
# Basic model structure: Y ~ Trt*Year + (1|Block) + (1|Sample_ID)
# Where Block = 1-8
# Where Sample_ID = Plot_ID (1-16)
# Used basic compound symmetry model (covariance structure)
#-------------------------------------------------------------------------------
source('analysis-functions.R')
#-------------------------------------------------------------------------------
# Libraries
library(data.table)
library(dplyr)
library(car)
library(emmeans)
library(ggplot2)
library(grid)
library(gridExtra)
library(lattice)
library(lmerTest)
library(lme4)
library(RColorBrewer)
library(rstudioapi)
library(viridis)
#-------------------------------------------------------------------------------
# Load directory paths
base_path    = dirname(getActiveDocumentContext()$path)
data_path    = paste(base_path, 'data', sep = '/')
figures_path = paste(base_path, 'figures', sep = '/')
tables_path  = paste(base_path, 'tables', sep = '/')
#-------------------------------------------------------------------------------
plfa_dt = fread(paste(data_path, 'AllYears_PLFA.csv', sep = '/'))
plfa_dt[, Year      := as.factor(Year)]
plfa_dt[, Block     := as.factor(Block)]
plfa_dt[, Sample_ID := as.factor(Sample_ID)]
str(plfa_dt)
#-------------------------------------------------------------------------------
# Total PLFA Biomass
variable = 'tot.biomass'
bwp.plfa = boxplot_check2(plfa_dt, 'tot.biomass')

plfa_dt[!is.na(tot.biomass), count := as.numeric(1L)]
plfa_dt[is.na(tot.biomass),  count := as.numeric(0L)]

summary_plfa_biomass = plfa_dt[, .(
  n                  = sum(count),
  mean        = round(mean(tot.biomass, na.rm = TRUE), digits = 1),
  sd                 = round(sd(tot.biomass, na.rm = TRUE), digits  = 1)),
  by                 = .(Trt, Year)]
summary_plfa_biomass[, SE := round(sqrt(sd/n), digits = 1)]

PLFA_m  = repeated_MM2(plfa_dt, variable)
#-------------------------------------------------------------------------------
# Subset data by key biomarkers

# Individual biomarkers
# 14:0 iso biomarker (gram-positive)
# 15:0 iso (gram-positive)
# 15:0 anteiso (gram-positive)
# 16:0 iso (gram-positive)
# 16:1 w7c (gram-negative)
# 16:1 w5c (AMF)
# 16:0 10-methyl biomarker (gram-positive, actinomycetes)
# 17:0 iso biomarker (gram-positive)
# 17:0 anteiso biomarker (gram-positive)
# 18:2 w6c biomarker (ectomycorrhizae and saprophytic fungi)
# 18:1 w9c biomarker (saprophytic fungi)
# 18:1 w7c biomarker (gram-negative bacteria)
# 18:0 10-methyl biomarker (gram-positive bacteria, actinomycetes)

# Grouped biomarkers

# Fungi: 18:2 w6c (ectomycorrhizae and saprophytic fungi) and 16:1 w5c (AMF)
# Bacteria total: 14:0 iso biomarker (gram-positive), 15:0 iso (gram-positive), 15:0 anteiso (gram-positive), 15:0,
#                 16:0 iso (gram-positive), 16:1 w7c (gram-negative), 17:0 iso biomarker (gram-positive)
#                 17:0 anteiso biomarker (gram-positive), 17:0, 17:1 w8c, 18:1 w5c, 18:1 w7c biomarker (gram-negative bacteria),
#                 16:0 10-methyl biomarker (gram-positive, actinomycetes), 17:0 10-methyl, 18:0 10-methyl biomarker (gram-positive bacteria, actinomycetes)
# Gram-positive: i14:0, i15:0, a15:0, i16:0, i17:0, a17:0, 10Me16:0, 10Me17:0, 10Me18:0
# Gram-negative: 16:1??7c, 17:1??8c, 18:1??7c, 18:1??5c
# Actinomycetes: 10Me16:0, 10Me17:0, 10Me18:0

plfa_grouped = plfa_dt[, .(
  bacteria   = sum(`019: 14:0 iso`, `031: 15:0 iso`,
                     `033: 15:0 anteiso`, `041: 15:0`,
                     `045: 16:0 iso`, `050: 16:1 w7c`,
                     `062: 17:0 anteiso`, `071: 17:0`,
                     `061: 17:0 iso`, `064: 17:1 w8c`,
                     `083: 18:1 w5c`, `081: 18:1 w7c`,
                     `057: 16:0 10-methyl`, `074: 17:0 10-methyl`,
                     `089: 18:0 10-methyl`, na.rm = TRUE),
  fungi      = sum(`078: 18:2 w6c`, `052: 16:1 w5c`, na.rm = TRUE),
  g_pos      = sum(`019: 14:0 iso`, `031: 15:0 iso`,
                   `033: 15:0 anteiso`, `045: 16:0 iso`,
                   `061: 17:0 iso`, `062: 17:0 anteiso`,
                   `057: 16:0 10-methyl`, `074: 17:0 10-methyl`,
                   `089: 18:0 10-methyl`, na.rm = TRUE),
  g_neg      = sum(`050: 16:1 w7c`, `064: 17:1 w8c`, 
                   `083: 18:1 w5c`, `081: 18:1 w7c`, na.rm = TRUE),
  actino     = sum(`057: 16:0 10-methyl`, `074: 17:0 10-methyl`,
                   `089: 18:0 10-methyl`, na.rm = TRUE)),
  by = .(Sample_ID, Trt, Block, Year)]
plfa_grouped[, fb_ratio   := (fungi/bacteria)]
head(plfa_grouped)
plfa_grouped[is.infinite(fb_ratio), fb_ratio := NA]
#-------------------------------------------------------------------------------
# Total bacteria
variable = 'bacteria'
bwp.bacteria = boxplot_check2(plfa_grouped, 'bacteria')

plfa_bacteria = plfa_grouped

plfa_bacteria[!is.na(bacteria), count := as.numeric(1L)]
plfa_bacteria[is.na(bacteria),  count := as.numeric(0L)]

summary_plfa_bacteria = plfa_bacteria[, .(
  n                  = sum(count),
  mean        = round(mean(bacteria, na.rm = TRUE), digits = 1),
  sd                 = round(sd(bacteria, na.rm = TRUE), digits  = 1)),
  by                 = .(Trt, Year)]
summary_plfa_bacteria[, SE := round(sqrt(sd/n), digits = 1)]

bacteria_m  = repeated_MM2(plfa_bacteria, variable)
# Total gram positive
variable = 'g_pos'
bwp.g_pos = boxplot_check2(plfa_grouped, 'g_pos')

plfa_g_pos = plfa_grouped

plfa_g_pos[!is.na(g_pos), count := as.numeric(1L)]
plfa_g_pos[is.na(g_pos),  count := as.numeric(0L)]

summary_plfa_g_pos = plfa_g_pos[, .(
  n                  = sum(count),
  mean        = round(mean(g_pos, na.rm = TRUE), digits = 1),
  sd                 = round(sd(g_pos, na.rm = TRUE), digits  = 1)),
  by                 = .(Trt, Year)]
summary_plfa_g_pos[, SE := round(sqrt(sd/n), digits = 1)]

g_pos_m  = repeated_MM2(plfa_g_pos, variable)
# Total gram negative
variable = 'g_neg'
bwp.g_neg = boxplot_check2(plfa_grouped, 'g_neg')

plfa_g_neg = plfa_grouped

plfa_g_neg[!is.na(g_neg), count := as.numeric(1L)]
plfa_g_neg[is.na(g_neg),  count := as.numeric(0L)]

summary_plfa_g_neg = plfa_g_neg[, .(
  n                  = sum(count),
  mean        = round(mean(g_neg, na.rm = TRUE), digits = 1),
  sd                 = round(sd(g_neg, na.rm = TRUE), digits  = 1)),
  by                 = .(Trt, Year)]
summary_plfa_g_neg[, SE := round(sqrt(sd/n), digits = 1)]

g_neg_m  = repeated_MM2(plfa_g_neg, variable)
# Total actino
variable = 'actino'
bwp.actino = boxplot_check2(plfa_grouped, 'actino')

plfa_actino = plfa_grouped

plfa_actino[!is.na(actino), count := as.numeric(1L)]
plfa_actino[is.na(actino),  count := as.numeric(0L)]

summary_plfa_actino = plfa_actino[, .(
  n                  = sum(count),
  mean        = round(mean(actino, na.rm = TRUE), digits = 1),
  sd                 = round(sd(actino, na.rm = TRUE), digits  = 1)),
  by                 = .(Trt, Year)]
summary_plfa_actino[, SE := round(sqrt(sd/n), digits = 1)]

actino_m  = repeated_MM2(plfa_actino, variable)
# Total fungi
variable = 'fungi'
bwp.fungi = boxplot_check2(plfa_grouped, 'fungi')

plfa_fungi = plfa_grouped

plfa_fungi[!is.na(fungi), count := as.numeric(1L)]
plfa_fungi[is.na(fungi),  count := as.numeric(0L)]

summary_plfa_fungi = plfa_fungi[, .(
  n                  = sum(count),
  mean        = round(mean(fungi, na.rm = TRUE), digits = 1),
  sd                 = round(sd(fungi, na.rm = TRUE), digits  = 1)),
  by                 = .(Trt, Year)]
summary_plfa_fungi[, SE := round(sqrt(sd/n), digits = 1)]

fungi_m  = repeated_MM2(plfa_fungi, variable)
# F:B
variable = 'fb_ratio'
bwp.fb_ratio = boxplot_check2(plfa_grouped, 'fb_ratio')

plfa_fb_ratio = plfa_grouped

plfa_fb_ratio[!is.na(fb_ratio), count := as.numeric(1L)]
plfa_fb_ratio[is.na(fb_ratio),  count := as.numeric(0L)]

summary_plfa_fb_ratio = plfa_fb_ratio[, .(
  n                  = sum(count),
  mean        = round(mean(fb_ratio, na.rm = TRUE), digits = 1),
  sd                 = round(sd(fb_ratio, na.rm = TRUE), digits  = 1)),
  by                 = .(Trt, Year)]
summary_plfa_fb_ratio[, SE := round(sqrt(sd/n), digits = 2)]

fb_ratio_m  = repeated_MM2(plfa_fb_ratio, variable)
#-------------------------------------------------------------------------------
# Output Tables
# Total PLFA biomass
sink(paste(tables_path, 'plfa_biomass_Repeated_Measures_Results.txt', sep = '/'))
print(PLFA_m)
sink()
sink(paste(tables_path, 'plfa_bacteria_Repeated_Measures_Results.txt', sep = '/'))
print(bacteria_m)
sink()
sink(paste(tables_path, 'plfa_g-pos_Repeated_Measures_Results.txt', sep = '/'))
print(g_pos_m)
sink()
sink(paste(tables_path, 'plfa_g-neg_Repeated_Measures_Results.txt', sep = '/'))
print(g_neg_m)
sink()
sink(paste(tables_path, 'plfa_actino_Repeated_Measures_Results.txt', sep = '/'))
print(actino_m)
sink()
sink(paste(tables_path, 'plfa_fungi_Repeated_Measures_Results.txt', sep = '/'))
print(fungi_m)
sink()
sink(paste(tables_path, 'plfa_fb-ratio_Repeated_Measures_Results.txt', sep = '/'))
print(fb_ratio_m)
sink()
#-------------------------------------------------------------------------------
f_labels = c('Post.2018'= '2018', 'Post.2019' = '2019', 'Post.2020' = '2020')
plfa_bio_gg = ggplot(plfa_dt[Year %in% c('Post.2018','Post.2019','Post.2020')],aes(x = Trt, y = tot.biomass, fill = Trt, group = Trt)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot(alpha = 0.9) +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  ylim(0.0,1560) +
  facet_grid(~ Year, labeller = labeller(Year = f_labels)) +
  scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.6, end = 0.2) +
  ylab(expression("PLFA Biomass (nmol g"^-1*")")) +
  xlab("Treatment") +
  theme_bw() +
  theme(text=element_text(size=7, color = 'black'),
        axis.title = element_text(size = 8, color = 'black'),
        strip.text.x = element_text(size = 7, color = 'black'),
        axis.text=element_text(size=7, color = 'black'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.title.y =  element_blank(),
        axis.title.x = element_blank(),
        legend.title=element_text(size = 7, color = 'black'),
        legend.position = "none",
        plot.title.position = "plot",  # This moves the title to align with plot edge
        plot.title = element_text(
          hjust = -0.005,  # Slight adjustment left of the plot
          vjust = -0.5,   # Slight adjustment above the plot
          size = 7       # Match your other text size if needed
        )) +
  guides(fill=guide_legend(title="Treatment")) +
  ggtitle('(a)')
plfa_bio_gg

# add P values
text = data.table(Trt = 'Compost', tot.biomass = 1500, lab = "n.s.",
                  Year  = factor(c("Post.2018"),levels = c("Post.2018","Post.2019","Post.2020")))
plfa_bio_gg = plfa_bio_gg + geom_text(data = text, label = "n.s.", 
                                  nudge_x = 0.5, size = 3)
text = data.table(Trt = 'Compost', tot.biomass = 1500, lab = "n.s.",
                  Year  = factor(c("Post.2019"),levels = c("Post.2018","Post.2019","Post.2020")))
plfa_bio_gg = plfa_bio_gg + geom_text(data = text, label = "n.s.", 
                                  nudge_x = 0.5, size = 3)
text = data.table(Trt = 'Compost', tot.biomass = 1500, lab = "P = 0.07",
                  Year  = factor(c("Post.2020"),levels = c("Post.2018","Post.2019","Post.2020")))
plfa_bio_gg = plfa_bio_gg + geom_text(data = text, label = "P = 0.07", 
                                  nudge_x = 0.5, size = 3)
plfa_bio_gg

bacteria_bio_gg = ggplot(plfa_grouped[Year %in% c('Post.2018','Post.2019','Post.2020')],aes(x = Trt, y = bacteria, fill = Trt, group = Trt)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot(alpha = 0.9) +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  ylim(0.0,850) +
  facet_grid(~ Year, labeller = labeller(Year = f_labels)) +
  scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.6, end = 0.2) +
  ylab(expression("PLFA Biomass (nmol g"^-1*")")) +
  xlab("Treatment") +
  theme_bw() +
  theme(text=element_text(size=7, color = 'black'),
        axis.title = element_text(size = 8, color = 'black'),
        strip.text.x = element_text(size = 7, color = 'black'),
        axis.text=element_text(size=7, color = 'black'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.title.y =  element_blank(),
        axis.title.x = element_blank(),
        legend.title=element_text(size = 7, color = 'black'),
        legend.position = "none",
        plot.title.position = "plot",  # This moves the title to align with plot edge
        plot.title = element_text(
          hjust = -0.005,  # Slight adjustment left of the plot
          vjust = -0.5,   # Slight adjustment above the plot
          size = 7       # Match your other text size if needed
        )) +
  guides(fill=guide_legend(title="Treatment")) +
  ggtitle('(b)')
bacteria_bio_gg

# add P values
text = data.table(Trt = 'Compost', bacteria = 800, lab = "n.s.",
                  Year  = factor(c("Post.2018"),levels = c("Post.2018","Post.2019","Post.2020")))
bacteria_bio_gg = bacteria_bio_gg + geom_text(data = text, label = "n.s.", 
                                      nudge_x = 0.5, size = 3)
text = data.table(Trt = 'Compost', bacteria = 800, lab = "n.s.",
                  Year  = factor(c("Post.2019"),levels = c("Post.2018","Post.2019","Post.2020")))
bacteria_bio_gg = bacteria_bio_gg + geom_text(data = text, label = "n.s.", 
                                      nudge_x = 0.5, size = 3)
text = data.table(Trt = 'Compost', bacteria = 800, lab = "P = 0.09",
                  Year  = factor(c("Post.2020"),levels = c("Post.2018","Post.2019","Post.2020")))
bacteria_bio_gg = bacteria_bio_gg + geom_text(data = text, label = "P = 0.09", 
                                      nudge_x = 0.5, size = 3)
bacteria_bio_gg

g_pos_bio_gg = ggplot(plfa_grouped[Year %in% c('Post.2018','Post.2019','Post.2020')],aes(x = Trt, y = g_pos, fill = Trt, group = Trt)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot(alpha = 0.9) +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  ylim(0.0,850) +
  facet_grid(~ Year, labeller = labeller(Year = f_labels)) +
  scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.6, end = 0.2) +
  ylab(expression("PLFA Biomass (nmol g"^-1*")")) +
  xlab("Treatment") +
  theme_bw() +
  theme(text=element_text(size=7, color = 'black'),
        axis.title = element_text(size = 8, color = 'black'),
        strip.text.x = element_text(size = 7, color = 'black'),
        axis.text=element_text(size=7, color = 'black'),
        axis.title.y =  element_blank(),
        axis.ticks.x = element_blank(),
        legend.title=element_text(size = 7, color = 'black'),
        legend.position = "none",
        plot.title.position = "plot",  # This moves the title to align with plot edge
        plot.title = element_text(
          hjust = -0.005,  # Slight adjustment left of the plot
          vjust = -0.5,   # Slight adjustment above the plot
          size = 7       # Match your other text size if needed
        )) +
  guides(fill=guide_legend(title="Treatment")) +
  ggtitle('(c)')
g_pos_bio_gg

# add P values
text = data.table(Trt = 'Compost', g_pos = 800, lab = "n.s.",
                  Year  = factor(c("Post.2018"),levels = c("Post.2018","Post.2019","Post.2020")))
g_pos_bio_gg = g_pos_bio_gg + geom_text(data = text, label = "n.s.", 
                                              nudge_x = 0.5, size = 3)
text = data.table(Trt = 'Compost', g_pos = 800, lab = "n.s.",
                  Year  = factor(c("Post.2019"),levels = c("Post.2018","Post.2019","Post.2020")))
g_pos_bio_gg = g_pos_bio_gg + geom_text(data = text, label = "n.s.", 
                                              nudge_x = 0.5, size = 3)
text = data.table(Trt = 'Compost', g_pos = 800, lab = "P = 0.06",
                  Year  = factor(c("Post.2020"),levels = c("Post.2018","Post.2019","Post.2020")))
g_pos_bio_gg = g_pos_bio_gg + geom_text(data = text, label = "P = 0.06", 
                                              nudge_x = 0.5, size = 3)
g_pos_bio_gg

actino_bio_gg = ggplot(plfa_grouped[Year %in% c('Post.2018','Post.2019','Post.2020')],aes(x = Trt, y = actino, fill = Trt, group = Trt)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot(alpha = 0.9) +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  ylim(0.0,850) +
  facet_grid(~ Year, labeller = labeller(Year = f_labels)) +
  scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.6, end = 0.2) +
  ylab(expression("PLFA Biomass (nmol g"^-1*")")) +
  xlab("Treatment") +
  theme_bw() +
  theme(text=element_text(size=7, color = 'black'),
        axis.title = element_text(size = 8, color = 'black'),
        strip.text.x = element_text(size = 7, color = 'black'),
        axis.text=element_text(size=7, color = 'black'),
        axis.title.y =  element_blank(),
        axis.ticks.x = element_blank(),
        legend.title=element_text(size = 7, color = 'black'),
        legend.position = "none",
        plot.title.position = "plot",  # This moves the title to align with plot edge
        plot.title = element_text(
          hjust = -0.005,  # Slight adjustment left of the plot
          vjust = -0.5,   # Slight adjustment above the plot
          size = 7       # Match your other text size if needed
        )) +
  guides(fill=guide_legend(title="Treatment")) +
  ggtitle('(d)')
actino_bio_gg

# add P values
text = data.table(Trt = 'Compost', actino = 800, lab = "n.s.",
                  Year  = factor(c("Post.2018"),levels = c("Post.2018","Post.2019","Post.2020")))
actino_bio_gg = actino_bio_gg + geom_text(data = text, label = "n.s.", 
                                        nudge_x = 0.5, size = 3)
text = data.table(Trt = 'Compost', actino = 800, lab = "n.s.",
                  Year  = factor(c("Post.2019"),levels = c("Post.2018","Post.2019","Post.2020")))
actino_bio_gg = actino_bio_gg + geom_text(data = text, label = "n.s.", 
                                        nudge_x = 0.5, size = 3)
text = data.table(Trt = 'Compost', actino = 800, lab = "P = 0.05",
                  Year  = factor(c("Post.2020"),levels = c("Post.2018","Post.2019","Post.2020")))
actino_bio_gg = actino_bio_gg + geom_text(data = text, label = "P = 0.05", 
                                        nudge_x = 0.5, size = 3)
actino_bio_gg

#-------------------------------------------------------------------------------
# Combine plots
plfa_p = grid.arrange(plfa_bio_gg, bacteria_bio_gg, g_pos_bio_gg, actino_bio_gg,
                         nrow = 2, ncol = 2)
ggsave("Figure3.tiff", plot = plfa_p, path = figures_path, width = 180, height = 180, units = "mm", dpi = 600)
