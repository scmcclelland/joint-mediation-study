# Plant Biomass and Community Analysis Script 
# Author: Shelby C McClelland
# Created:     20 December 2020
# Last Update: 02 August 2024
# Description: This file analyzes soil chemical property data.
#-------------------------------------------------------------------------------
# Analysis notes:
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
library(ggtext)
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
# Read in files
shoot_dt = fread(paste(data_path, 'AllYears_Shoot_Biomass.csv', sep = '/'))
shoot_dt[, Year      := as.factor(Year)]
shoot_dt[, Block     := as.factor(Block)]
shoot_dt[, Sample_ID := as.factor(Sample_ID)]
str(shoot_dt)
root_dt  = fread(paste(data_path, 'AllYears_Root_Biomass.csv', sep = '/'))
root_dt[, Year      := as.factor(Year)]
root_dt[, Block     := as.factor(Block)]
root_dt[, Sample_ID := as.factor(Sample_ID)]
# Units: Mg C ha-1
div_dt  = fread(paste(data_path, 'AllYears_Plant_Diversity.csv', sep = '/'))
div_dt[, Year      := as.factor(Year)]
div_dt[, Block     := as.factor(Block)]
div_dt[, Sample_ID := as.factor(Sample_ID)]
relabun_dt  = fread(paste(data_path, 'AllYears_Plant_lnRR.csv', sep = '/'))
relabun_dt[, Year      := as.factor(Year)]
relabun_dt[, Block     := as.factor(Block)]
#-------------------------------------------------------------------------------
# Roots (total biomass)
variable = 'total'
bwp.totalroot = boxplot_check(root_dt, 'total', 'Trt')

totalroot_dt  = root_dt[, c('Sample_ID', 'Trt', 'Block', 'Depth', 'Year', 'total')]
totalroot_dt[!is.na(total), count := as.numeric(1L)]
totalroot_dt[is.na(total),  count := as.numeric(0L)]

summary_totalroot = totalroot_dt[, .(
  n        = sum(count),
  mean_total  = round(mean(total, na.rm = TRUE), digits = 1),
  sd       = round(sd(total, na.rm = TRUE), digits  = 1)),
  by       = .(Trt, Year, Depth)]
summary_totalroot[, SE := round(sqrt(sd/n), digits = 1)]

totalroot_0_10cm = repeated_MM(totalroot_dt, variable, '0to10')
totalroot_10_20cm = repeated_MM(totalroot_dt, variable, '10to20')

# Roots (coarse biomass)
variable = 'coarse'
bwp.coarseroot = boxplot_check(root_dt, 'coarse', 'Trt')

coarseroot_dt  = root_dt[, c('Sample_ID', 'Trt', 'Block', 'Depth', 'Year', 'coarse')]
coarseroot_dt[!is.na(coarse), count := as.numeric(1L)]
coarseroot_dt[is.na(coarse),  count := as.numeric(0L)]

summary_coarseroot = coarseroot_dt[, .(
  n        = sum(count),
  mean_coarse  = round(mean(coarse, na.rm = TRUE), digits = 1),
  sd       = round(sd(coarse, na.rm = TRUE), digits  = 1)),
  by       = .(Trt, Year, Depth)]
summary_coarseroot[, SE := round(sqrt(sd/n), digits = 1)]

coarseroot_0_10cm = repeated_MM(coarseroot_dt, variable, '0to10')
coarseroot_10_20cm = repeated_MM(coarseroot_dt, variable, '10to20')

# Roots (fine biomass)
variable = 'fine'
bwp.fineroot = boxplot_check(root_dt, 'fine', 'Trt')

fineroot_dt  = root_dt[, c('Sample_ID', 'Trt', 'Block', 'Depth', 'Year', 'fine')]
fineroot_dt[!is.na(fine), count := as.numeric(1L)]
fineroot_dt[is.na(fine),  count := as.numeric(0L)]

summary_fineroot = fineroot_dt[, .(
  n        = sum(count),
  mean_fine  = round(mean(fine, na.rm = TRUE), digits = 1),
  sd       = round(sd(fine, na.rm = TRUE), digits  = 1)),
  by       = .(Trt, Year, Depth)]
summary_fineroot[, SE := round(sqrt(sd/n), digits = 1)]

fineroot_0_10cm = repeated_MM(fineroot_dt, variable, '0to10')
fineroot_10_20cm = repeated_MM(fineroot_dt, variable, '10to20')
#-------------------------------------------------------------------------------
# Output Tables
# Total root biomass
sink(paste(tables_path, 'total_root_biomass_Repeated_Measures_Results.txt', sep = '/'))
print(totalroot_0_10cm)
print(totalroot_10_20cm)
sink()
# Coarse root biomass
sink(paste(tables_path, 'coarse_root_biomass_Repeated_Measures_Results.txt', sep = '/'))
print(coarseroot_0_10cm)
print(coarseroot_10_20cm)
sink()
# Fine root biomass
sink(paste(tables_path, 'fine_root_biomass_Repeated_Measures_Results.txt', sep = '/'))
print(fineroot_0_10cm)
print(fineroot_10_20cm)
sink()
#-------------------------------------------------------------------------------
# Create Root Biomass Figures
f_labels = c('2018.su'= '2018', '2019' = '2019', '2020.su' = '2020')
t_root_bio_gg <- ggplot(root_dt[Year %in% c('2018.su','2019','2020.su')],aes(x = Trt, y = total, fill = Trt, group = Trt)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot(alpha = 0.9) +
  stat_summary(fun=mean, geom="point", shape=23, size=4) +
  ylim(0.0,4.5) +
  facet_grid(~ Year, labeller = labeller(Year = f_labels)) +
  scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.6, end = 0.2) +
  ylab(expression("Biomass Carbon (Mg C ha"^-1*")")) +
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=13, color = 'black'),
        strip.text.x = element_text(size = 13, color = 'black'),
        axis.text=element_text(size=13, color = 'black'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y =  element_blank(),
        axis.title.x = element_blank(),
        legend.title=element_text(size = 13, color = 'black'),
        legend.position = "none") +
  guides(fill=guide_legend(title="Treatment"))
t_root_bio_gg

c_root_bio_gg <- ggplot(root_dt[Year %in% c('2018.su','2019','2020.su')],aes(x = Trt, y = coarse, fill = Trt, group = Trt)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot(alpha = 0.9) +
  stat_summary(fun=mean, geom="point", shape=23, size=4) +
  ylim(0.0,4.5) +
  facet_grid(~ Year, labeller = labeller(Year = f_labels)) +
  scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.6, end = 0.2) +
  ylab(expression("Biomass Carbon (Mg C ha"^-1*")")) +
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=13, color = 'black'),
        strip.text.x = element_text(size = 13, color = 'black'),
        axis.text    = element_text(size=13, color = 'black'),
        axis.title.y =  element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 13, color = 'black'),
        legend.position = "none") +
  guides(fill=guide_legend(title="Treatment"))
c_root_bio_gg

f_root_bio_gg <- ggplot(root_dt[Year %in% c('2018.su','2019','2020.su')],aes(x = Trt, y = fine, fill = Trt, group = Trt)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot(alpha = 0.9) +
  stat_summary(fun=mean, geom="point", shape=23, size=4) +
  ylim(0.0,4.5) +
  facet_grid(~ Year, labeller = labeller(Year = f_labels)) +
  scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.6, end = 0.2) +
  ylab(expression("Biomass Carbon (Mg C ha"^-1*")")) +
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=13, color = 'black'),
        strip.text.x = element_text(size = 13, color = 'black'),
        axis.text=element_text(size=13, color = 'black'),
        axis.title.y =  element_blank(),
        axis.ticks.x = element_blank(),
        legend.title=element_text(size = 13, color = 'black'),
        legend.position = "none") +
  guides(fill=guide_legend(title="Treatment"))
f_root_bio_gg
#-------------------------------------------------------------------------------
# Shoot (biomass C)
variable = 'Biomass_C'
bwp.shoot = boxplot_check2(shoot_dt, 'Biomass_C')

shoot_dt[!is.na(Biomass_C), count := as.numeric(1L)]
shoot_dt[is.na(Biomass_C),  count := as.numeric(0L)]

summary_shoot = shoot_dt[, .(
  n           = sum(count),
  mean_shootC = round(mean(Biomass_C, na.rm = TRUE), digits = 1),
  sd          = round(sd(Biomass_C, na.rm = TRUE), digits  = 1)),
  by          = .(Trt, Year)]
summary_shoot[, SE := round(sqrt(sd/n), digits = 1)]

BiomassC_m  = repeated_MM2(shoot_dt, variable)
#-------------------------------------------------------------------------------
# Output Tables
# Shoot biomass
sink(paste(tables_path, 'shoot_biomass_Repeated_Measures_Results.txt', sep = '/'))
print(BiomassC_m)
sink()
#-------------------------------------------------------------------------------
# Create Shoot Biomass Figures
shoot_bio_gg <- ggplot(shoot_dt,aes(x = Trt, y = Biomass_C, fill = Trt, group = Trt)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot(alpha = 0.9) +
  stat_summary(fun=mean, geom="point", shape=23, size=4) +
  ylim(0.0,4.5) +
  facet_grid(~ Year, labeller = label_parsed) +
  scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.6, end = 0.2) +
  ylab(expression("Biomass Carbon (Mg C ha"^-1*")")) +
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=13, color = 'black'),
        strip.text.x = element_text(size = 13, color = 'black'),
        axis.text=element_text(size=13, color = 'black'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y =  element_blank(),
        axis.title.x = element_blank(),
        legend.title=element_text(size = 13, color = 'black'),
        legend.position = "none") +
  guides(fill=guide_legend(title="Treatment"))
shoot_bio_gg
#-------------------------------------------------------------------------------
# Combine plots
yleft = textGrob(expression(paste("Biomass Carbon (Mg C ha"^-1*")")), 
                  rot = 90, gp = gpar(fontsize = 16))
biomass_p = grid.arrange(shoot_bio_gg, t_root_bio_gg, c_root_bio_gg, f_root_bio_gg,
                         nrow = 2, ncol = 2,
                         left = yleft)
ggsave("EcolLetters_Biomass.tiff", plot = biomass_p, path = figures_path, width = 11, height = 8, units = "in", dpi = 300)
#-------------------------------------------------------------------------------
# Plant Community Analyses
# Shannon H
variable = 'Avg_H'
bwp.H = boxplot_check2(div_dt, 'Avg_H')

div_H_dt  = div_dt[, c('Sample_ID', 'Trt', 'Block', 'Year', 'Avg_H')]
div_H_dt[!is.na(Avg_H), count := as.numeric(1L)]
div_H_dt[is.na(Avg_H),  count := as.numeric(0L)]

summary_div_H = div_H_dt[, .(
  n           = sum(count),
  mean_total  = round(mean(Avg_H, na.rm = TRUE), digits = 1),
  sd          = round(sd(Avg_H, na.rm = TRUE), digits  = 1)),
  by          = .(Trt, Year)]
summary_div_H[, SE := round(sqrt(sd/n), digits = 1)]

div_H_MM   = repeated_MM2(div_H_dt, variable)
# Evenness
variable = 'Avg_Evenness'
bwp.Ev   = boxplot_check2(div_dt, 'Avg_Evenness')

div_Ev_dt  = div_dt[, c('Sample_ID', 'Trt', 'Block', 'Year', 'Avg_Evenness')]
div_Ev_dt[!is.na(Avg_Evenness), count := as.numeric(1L)]
div_Ev_dt[is.na(Avg_Evenness),  count := as.numeric(0L)]

summary_div_Ev = div_Ev_dt[, .(
  n            = sum(count),
  mean_total   = round(mean(Avg_Evenness, na.rm = TRUE), digits = 1),
  sd           = round(sd(Avg_Evenness, na.rm = TRUE), digits  = 1)),
  by           = .(Trt, Year)]
summary_div_Ev[, SE := round(sqrt(sd/n), digits = 1)]

div_Ev_MM   = repeated_MM2(div_Ev_dt, variable)
# Richness
variable = 'Avg_Richness'
bwp.Ri = boxplot_check2(div_dt, 'Avg_Richness')

div_Ri_dt  = div_dt[, c('Sample_ID', 'Trt', 'Block', 'Year', 'Avg_Richness')]
div_Ri_dt[!is.na(Avg_Richness), count := as.numeric(1L)]
div_Ri_dt[is.na(Avg_Richness),  count := as.numeric(0L)]

summary_div_Ri = div_Ri_dt[, .(
  n        = sum(count),
  mean_total  = round(mean(Avg_Richness, na.rm = TRUE), digits = 1),
  sd       = round(sd(Avg_Richness, na.rm = TRUE), digits  = 1)),
  by       = .(Trt, Year)]
summary_div_Ri[, SE := round(sqrt(sd/n), digits = 1)]

div_Ri_MM   = repeated_MM2(div_Ri_dt, variable)
#-------------------------------------------------------------------------------
# Output Tables
# Plant diversity metrics
sink(paste(tables_path, 'Shannon_Plant_Diversity_Repeated_Measures_Results.txt', sep = '/'))
print(div_H_MM)
sink()
sink(paste(tables_path, 'Plant_Diversity_Evenness_Repeated_Measures_Results.txt', sep = '/'))
print(div_Ev_MM)
sink()
sink(paste(tables_path, 'Plant_Diversity_Richness_Repeated_Measures_Results.txt', sep = '/'))
print(div_Ri_MM)
sink()
#-------------------------------------------------------------------------------
# Create Plant Diversity Figures
div_H_gg <- ggplot(div_dt,aes(x = Trt, y = Avg_H, fill = Trt, group = Trt)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot(alpha = 0.9) +
  stat_summary(fun=mean, geom="point", shape=23, size=4) +
  ylim(0.0,1.5) +
  facet_grid(~ Year) +
  scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.6, end = 0.2) +
  ylab("Shannon Diversity (H)") +
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=13, color = 'black'),
        strip.text.x = element_text(size = 13, color = 'black'),
        axis.text=element_text(size=13, color = 'black'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.title.y =  element_blank(),
        axis.title.x = element_blank(),
        legend.title=element_text(size = 13, color = 'black'),
        legend.position = "none") +
  guides(fill=guide_legend(title="Treatment"))
div_H_gg

div_Ev_gg <- ggplot(div_dt,aes(x = Trt, y = Avg_Evenness, fill = Trt, group = Trt)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot(alpha = 0.9) +
  stat_summary(fun=mean, geom="point", shape=23, size=4) +
  ylim(0.0,1.5) +
  facet_grid(~ Year) +
  scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.6, end = 0.2) +
  ylab("Plant Species Evenness") +
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=13, color = 'black'),
        strip.text.x = element_text(size = 13, color = 'black'),
        axis.text=element_text(size=13, color = 'black'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.title.y =  element_blank(),
        axis.title.x = element_blank(),
        legend.title=element_text(size = 13, color = 'black'),
        legend.position = "none") +
  guides(fill=guide_legend(title="Treatment"))
div_Ev_gg

div_Ri_gg <- ggplot(div_dt,aes(x = Trt, y = Avg_Richness, fill = Trt, group = Trt)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot(alpha = 0.9) +
  stat_summary(fun=mean, geom="point", shape=23, size=4) +
  ylim(0.0,5) +
  facet_grid(~ Year) +
  scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.6, end = 0.2) +
  ylab("Plant Species Richness") +
  xlab("Treatment") +
  theme_classic() +
  theme(text=element_text(size=13, color = 'black'),
        strip.text.x = element_text(size = 13, color = 'black'),
        axis.text=element_text(size=13, color = 'black'),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.title.y =  element_blank(),
        axis.title.x = element_blank(),
        legend.title=element_text(size = 13, color = 'black'),
        legend.position = "none") +
  guides(fill=guide_legend(title="Treatment"))
div_Ri_gg
#-------------------------------------------------------------------------------
# Combine plots
div_p = grid.arrange(div_H_gg, div_Ev_gg, div_Ri_gg,
                         nrow = 3, ncol = 1)
ggsave("EcolLetters_PlantDiversity.tiff", plot = div_p, path = figures_path, width = 6, height = 9, units = "in", dpi = 300)
#-------------------------------------------------------------------------------
# Relative Abundance 
# aggregate inermis and commutatus to spp
relabun_dt[Species %like% 'Bromus', Species := 'Bromus spp']
relabun_dt[Species %in% 'Bromus spp', Rel.Abund := sum(Rel.Abund), by = .(Trt, Block, Year, Species)]
relabun_dt = unique(relabun_dt)

# mean response by year
relabund_m_dt   = relabun_dt[, .(
  Rel.Abund  = round(mean(Rel.Abund, na.rm = TRUE), digits = 1)),
  by       = .(Trt, Block, Species)]

# Transform from long to wide format
relabund_w_dt = dcast(relabund_m_dt, Species+Block ~ Trt,
                     value.var = 'Rel.Abund')
# Compute log response
relabund_w_dt[, lnRR := log(Compost/Control)]
relabund_w_dt = relabund_w_dt[complete.cases(relabund_w_dt)]
relabund_w_dt = relabund_w_dt[!lnRR == Inf & !lnRR == -Inf,]
relabund_w_dt[, count := as.numeric(1L)]
relabund_w_dt[, count := sum(count), by = Species]
relabund_w_dt = relabund_w_dt[count > 1,]

# one sample t-tests
bro.spp = t.test(relabund_w_dt[Species %in% 'Bromus spp', lnRR])    # p = 0.08
con.arv = t.test(relabund_w_dt[Species %like% 'Convolvulus', lnRR]) # NS
dac.glo = t.test(relabund_w_dt[Species %like% 'Dactylis', lnRR])    # p = 0.07
fes.aru = t.test(relabund_w_dt[Species %like% 'Festuca', lnRR])     # NS
med.sat = t.test(relabund_w_dt[Species %like% 'Medicago', lnRR])    # p = 0.10

# create table of t-test results, percent change
lnRR_dt = data.table(Species = unique(relabund_w_dt[, Species]), 
                     m_lnRR  = c(bro.spp$estimate, con.arv$estimate, dac.glo$estimate, fes.aru$estimate, med.sat$estimate),
                     l_CI    = c(bro.spp$conf.int[1], con.arv$conf.int[1], dac.glo$conf.int[1], fes.aru$conf.int[1], med.sat$conf.int[1]),
                     u_CI    = c(bro.spp$conf.int[2], con.arv$conf.int[2], dac.glo$conf.int[2], fes.aru$conf.int[2], med.sat$conf.int[2]))
perc_dt = copy(lnRR_dt)
perc_dt[, m_lnRR := (exp(m_lnRR)-1)*100]
perc_dt[, l_CI   := (exp(l_CI)-1)*100]
perc_dt[, u_CI   := (exp(u_CI)-1)*100]
#-------------------------------------------------------------------------------
# Output Tables
# Plant diversity metrics
sink(paste(tables_path, 'lnRR_Plant_Relative_Abundances.txt', sep = '/'))
print(bro.spp)
print(con.arv)
print(dac.glo)
print(fes.aru)
print(med.sat)
sink()
#-------------------------------------------------------------------------------
# Plot Log Response Figure
# easier to see with Festuca removed
lnRR_gg = ggplot(lnRR_dt[!Species %like% 'Festuca'],aes(x = Species, y = m_lnRR)) +
  geom_point(size = 3, color = 'black') + 
  geom_errorbar(aes(ymin = l_CI, ymax = u_CI), width = .2,
                position=position_dodge(0.05)) +
  coord_flip() + 
  geom_hline(yintercept = 0) +
  # scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.6, end = 0.2) +
  ylab("lnRR") +
  theme_bw() +
  theme(text=element_text(size=13, color = 'black'),
        strip.text.x = element_text(size = 13, color = 'black'), 
        axis.text.y = element_text(face = "italic"),
        axis.title.y = element_blank(),
        axis.text=element_text(size=13, color = 'black'),
        legend.position = "none")
lnRR_gg
ggsave("EcolLetters_lnRR.tiff", plot = lnRR_gg, path = figures_path, width = 5, height = 5, units = "in", dpi = 300)
