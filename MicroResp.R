# MicroResp Analysis Script 
# Author: Shelby C McClelland
# Created:     20 December 2020
# Last Update: 22 February 2025
# Description: This file analyzes MicroResp data.
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
mr_dt = fread(paste(data_path, 'AllYears_MicroResp.csv', sep = '/'))
mr_dt[, Year      := as.factor(Year)]
mr_dt[, Block     := as.factor(Block)]
mr_dt[, Sample_ID := as.factor(Sample_ID)]
str(mr_dt)

mre_dt = fread(paste(data_path, 'AllYears_MREvenness.csv', sep = '/'))
mre_dt[, Year      := as.factor(Year)]
mre_dt[, Block     := as.factor(Block)]
mre_dt[, Sample_ID := as.factor(Sample_ID)]
str(mre_dt)
#-------------------------------------------------------------------------------
# Substrate Analysis
variable = 'CO2'
bwp.mr.glu  = boxplot_check2(mr_dt[Substrate %in% 'Glucose'], 'CO2')
bwp.mr.cell = boxplot_check2(mr_dt[Substrate %in% 'Cellulose'], 'CO2')
bwp.mr.xyl  = boxplot_check2(mr_dt[Substrate %in% 'Xylose'], 'CO2')
bwp.mr.gluc = boxplot_check2(mr_dt[Substrate %in% 'Glucosamine'], 'CO2')
bwp.mr.lig  = boxplot_check2(mr_dt[Substrate %in% 'Lignin'], 'CO2')
bwp.mr.h2o  = boxplot_check2(mr_dt[Substrate %in% 'Water'], 'CO2')

mr_dt[!is.na(CO2), count := as.numeric(1L)]
mr_dt[is.na(CO2),  count := as.numeric(0L)]

summary_mr = mr_dt[, .(
  n                  = sum(count),
  mean_CO2           = round(mean(CO2, na.rm = TRUE), digits = 1),
  sd                 = round(sd(CO2, na.rm = TRUE), digits  = 1)),
  by                 = .(Trt, Year, Substrate)]
summary_mr[, SE := round(sqrt(sd/n), digits = 1)]

# Log transform data to meet normality assumptions
mr_dt[, CO2_log := log(CO2)]
variable = 'CO2_log'

MR_m      = repeated_MM3(mr_dt, variable)
#-------------------------------------------------------------------------------
# Evenness
variable    = 'MR.Func.Ev'
bwp.mr.ev  = boxplot_check2(mre_dt, 'MR.Func.Ev')

mre_dt[!is.na(MR.Func.Ev), count := as.numeric(1L)]
mre_dt[is.na(MR.Func.Ev),  count := as.numeric(0L)]

summary_mre  = mre_dt[, .(
  n                  = sum(count),
  mean_Ev            = round(mean(MR.Func.Ev, na.rm = TRUE), digits = 1),
  sd                 = round(sd(MR.Func.Ev, na.rm = TRUE), digits  = 1)),
  by                 = .(Trt, Year)]
summary_mre[, SE := round(sqrt(sd/n), digits = 1)]

MR_ev_m      = repeated_MM2(mre_dt, variable)
#-------------------------------------------------------------------------------
# Total Respiration
variable    = 'Tot.Resp'
bwp.mr.tr  = boxplot_check2(mre_dt, 'Tot.Resp')

mre_dt[!is.na(Tot.Resp), count := as.numeric(1L)]
mre_dt[is.na(Tot.Resp),  count := as.numeric(0L)]

summary_mre  = mre_dt[, .(
  n                  = sum(count),
  mean_TR            = round(mean(Tot.Resp, na.rm = TRUE), digits = 1),
  sd                 = round(sd(Tot.Resp, na.rm = TRUE), digits  = 1)),
  by                 = .(Trt, Year)]
summary_mre[, SE := round(sqrt(sd/n), digits = 1)]

MR_tr_m      = repeated_MM2(mre_dt, variable)
#-------------------------------------------------------------------------------
# Output Tables
# Substrate Analysis
sink(paste(tables_path, 'MR_substrates_Repeated_Measures_Results.txt', sep = '/'))
print(MR_m)
sink()
# Evenness Analysis
sink(paste(tables_path, 'MR_evenness_Repeated_Measures_Results.txt', sep = '/'))
print(MR_ev_m)
sink()
# Total Respiration
sink(paste(tables_path, 'MR_totalrespiration_Repeated_Measures_Results.txt', sep = '/'))
print(MR_tr_m)
sink()
#-------------------------------------------------------------------------------
# PLOT
# set order
mr_dt$Substrate <- factor(mr_dt$Substrate, levels = c("Lignin", "Glucosamine", "Xylose", "Cellulose", "Glucose", "Water"))
mr_dt$Substrate

MR_substr_gg = ggplot(mr_dt[!Substrate %in% 'Water'],aes(x = Substrate, y = CO2, 
                                group = interaction(Trt, Substrate),
                                 fill = Trt)) +
  stat_boxplot(geom = "errorbar") +  
  geom_boxplot(alpha = 0.9) +
  stat_summary(fun=mean, geom="point", shape=23, size=2, 
               position = position_dodge(width = 0.75), 
               aes(group = interaction(Trt, Substrate))) +
  facet_wrap(~ Year, labeller = label_parsed, nrow = 3, ncol = 1) +
  scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.6, end = 0.2) +
  ylab(expression("Substrate Induced Respiration ("*mu*"g CO"[2]*"-C g"^-1*" h"^-1*")")) +
  xlab("Substrate") +
  theme_bw() +
  theme(text=element_text(size=7, color = 'black'),
        axis.title = element_text(size = 8, color = 'black'),
        strip.text.x = element_text(size = 7, color = 'black'),
        axis.text    = element_text(size=7, color = 'black'),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 7, color = 'black'),
        legend.position = "none") +
  guides(fill=guide_legend(title="Treatment"))
MR_substr_gg

# add P-values by substrate x year
# glucose
text = data.table(
  Year = unique(mr_dt$Year),  # year levels
  Substrate = "Glucose",      # substrate position for the label
  CO2   = c(5.5, 5, 0.5),       # position
  Trt   = 'Compost',          # treatment
  label = c('P = 0.06', 'n.s.', 'n.s.')
)
MR_substr_gg = MR_substr_gg +
  geom_text(data = text,
            aes(x = Substrate, y = CO2, label = label, group = NULL),
            inherit.aes = FALSE,
            size = 3    # Adjust size to match your theme
  ) 
# cellulose
text = data.table(
  Year = unique(mr_dt$Year),  # year levels
  Substrate = "Cellulose",      # substrate position for the label
  CO2   = c(5.5, 5, 5),       # position
  Trt   = 'Compost',          # treatment
  label = c('n.s.', 'n.s.', 'n.s.')
)
MR_substr_gg = MR_substr_gg +
  geom_text(data = text,
            aes(x = Substrate, y = CO2, label = label, group = NULL),
            inherit.aes = FALSE,
            size = 3    # Adjust size to match your theme
  ) 
# xylose
text = data.table(
  Year = unique(mr_dt$Year),  # year levels
  Substrate = "Xylose",      # substrate position for the label
  CO2   = c(5.5, 5, 5),       # position
  Trt   = 'Compost',          # treatment
  label = c('n.s.', 'n.s.', 'P = 0.06')
)
MR_substr_gg = MR_substr_gg +
  geom_text(data = text,
            aes(x = Substrate, y = CO2, label = label, group = NULL),
            inherit.aes = FALSE,
            size = 3    # Adjust size to match your theme
  ) 
# glucosamine
text = data.table(
  Year = unique(mr_dt$Year),  # year levels
  Substrate = "Glucosamine",      # substrate position for the label
  CO2   = c(5.5, 5, 5),       # position
  Trt   = 'Compost',          # treatment
  label = c('n.s.', 'n.s.', 'n.s.')
)
MR_substr_gg = MR_substr_gg +
  geom_text(data = text,
            aes(x = Substrate, y = CO2, label = label, group = NULL),
            inherit.aes = FALSE,
            size = 3    # Adjust size to match your theme
  ) 
# lignin
text = data.table(
  Year = unique(mr_dt$Year),  # year levels
  Substrate = "Lignin",      # substrate position for the label
  CO2   = c(5.5, 5, 5),       # position
  Trt   = 'Compost',          # treatment
  label = c('n.s.', 'n.s.', 'n.s.')
)
MR_substr_gg = MR_substr_gg +
  geom_text(data = text,
            aes(x = Substrate, y = CO2, label = label, group = NULL),
            inherit.aes = FALSE,
            size = 3    # Adjust size to match your theme
  ) 
MR_substr_gg
ggsave("Figure4.tiff", plot = MR_substr_gg, path = figures_path, width = 88, height = 180, units = "mm", dpi = 600)
