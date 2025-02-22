# Soil Chemical Properties Analysis Script 
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
soilchem_dt = fread(paste(data_path, 'AllYears_SoilNutrients.csv', sep = '/'))
str(soilchem_dt)
soilchem_dt[, Block     := as.factor(Block)]
soilchem_dt[, Sample_ID := as.factor(Sample_ID)]
str(soilchem_dt)
# Units:
# Bulk density g cm-3
# Potassium: ppm
# Mehlich P: ppm

soilC_dt = fread(paste(data_path, 'AllYears_SOC_ESM.csv', sep = '/'))
# Units
# Mg C ha-1
#-------------------------------------------------------------------------------
#pH
bwp.pH = boxplot_check(soilchem_dt, 'pH', 'Trt')

soilpH_dt  = soilchem_dt[, c('Sample_ID', 'Trt', 'Block', 'Depth', 'Year', 'pH')]
soilpH_dt[!is.na(pH), count := as.numeric(1L)]
soilpH_dt[is.na(pH),  count := as.numeric(0L)]

summary_pH = soilpH_dt[, .(
  n        = sum(count),
  mean_pH  = round(mean(pH, na.rm = TRUE), digits = 1),
  sd       = round(sd(pH, na.rm = TRUE), digits  = 1)),
  by       = .(Trt, Year, Depth)]
summary_pH[, SE := round(sqrt(sd/n), digits = 1)]

variable = 'pH'
pH_0_10cm  = repeated_MM(soilpH_dt, variable, '0to10')
pH_10_20cm = repeated_MM(soilpH_dt, variable, '10to20')
pH_20_50cm = repeated_MM(soilpH_dt, variable, '20to50')
#-------------------------------------------------------------------------------
# Mehlich P 
bwp.mehlichP = boxplot_check(soilchem_dt, 'Mehlich.P.ppm', 'Trt', 'Depth', 'Year')
bwp.mehlichP

soilmehlichP_dt  = soilchem_dt[, c('Sample_ID', 'Trt', 'Block', 'Depth', 'Year', 'Mehlich.P.ppm')]
soilmehlichP_dt[!is.na(Mehlich.P.ppm), count := as.numeric(1L)]
soilmehlichP_dt[is.na(Mehlich.P.ppm),  count := as.numeric(0L)]

summary_mehlichP = soilmehlichP_dt[, .(
  n              = sum(count),
  mean_mehlichP  = round(mean(Mehlich.P.ppm, na.rm = TRUE), digits = 1),
  sd             = round(sd(Mehlich.P.ppm, na.rm = TRUE), digits  = 1)),
  by             = .(Trt, Year, Depth)]
summary_mehlichP[, SE := round(sqrt(sd/n), digits = 1)]

variable = 'Mehlich.P.ppm'
mehlichP_0_10cm  = repeated_MM(soilmehlichP_dt, variable, '0to10')
mehlichP_10_20cm = repeated_MM(soilmehlichP_dt, variable, '10to20')
mehlichP_20_50cm = MM(soilmehlichP_dt, variable, '20to50')
#-------------------------------------------------------------------------------
# Potassium 
bwp.Potassium = boxplot_check(soilchem_dt, 'Potassium.ppm', 'Trt', 'Depth', 'Year')
bwp.Potassium

soilPotassium_dt  = soilchem_dt[, c('Sample_ID', 'Trt', 'Block','Depth', 'Year', 'Potassium.ppm')]
soilPotassium_dt[!is.na(Potassium.ppm), count := as.numeric(1L)]
soilPotassium_dt[is.na(Potassium.ppm),  count := as.numeric(0L)]

summary_Potassium = soilPotassium_dt[, .(
  n              = sum(count),
  mean_Potassium  = round(mean(Potassium.ppm, na.rm = TRUE), digits = 1),
  sd             = round(sd(Potassium.ppm, na.rm = TRUE), digits  = 1)),
  by             = .(Trt, Year, Depth)]
summary_Potassium[, SE := round(sqrt(sd/n), digits = 1)]

variable = 'Potassium.ppm'
Potassium_0_10cm  = repeated_MM(soilPotassium_dt, variable, '0to10')
Potassium_10_20cm = repeated_MM(soilPotassium_dt, variable, '10to20')
Potassium_20_50cm = MM(soilPotassium_dt, variable, '20to50')
#-------------------------------------------------------------------------------
#Soil organic carbon
bwp.SOC = boxplot_check(soilC_dt, 'soc.Mg.ha', 'Trt')

soilSOC_dt  = soilC_dt[, c('Sample_ID', 'Trt', 'Block', 'Depth', 'Year', 'soc.Mg.ha')]
soilSOC_dt[!is.na(soc.Mg.ha), count := as.numeric(1L)]
soilSOC_dt[is.na(soc.Mg.ha),  count := as.numeric(0L)]

summary_SOC = soilSOC_dt[, .(
  n        = sum(count),
  mean_SOC  = round(mean(soc.Mg.ha, na.rm = TRUE), digits = 1),
  sd       = round(sd(soc.Mg.ha, na.rm = TRUE), digits  = 1)),
  by       = .(Trt, Year, Depth)]
summary_SOC[, SE := round(sqrt(sd/n), digits = 1)]

variable = 'soc.Mg.ha'
SOC_0_10cm = repeated_MM(soilSOC_dt, variable, '0to10')
SOC_10_20cm = repeated_MM(soilSOC_dt, variable, '10to20')
#-------------------------------------------------------------------------------
# Output Tables
#pH
sink(paste(tables_path, 'pH_Repeated_Measures_Results.txt', sep = '/'))
print(pH_0_10cm)
print(pH_10_20cm)
print(pH_20_50cm)
sink()
# Mehlich P
sink(paste(tables_path, 'Phosphorus_Repeated_Measures_Results.txt', sep = '/'))
print(mehlichP_0_10cm)
print(mehlichP_10_20cm)
print(mehlichP_20_50cm)
sink()
# Potassium
sink(paste(tables_path, 'Potassium_Repeated_Measures_Results.txt', sep = '/'))
print(Potassium_0_10cm)
print(Potassium_10_20cm)
print(Potassium_20_50cm)
sink()
# Soil Organic Carbon
sink(paste(tables_path, 'SOC_Repeated_Measures_Results.txt', sep = '/'))
print(SOC_0_10cm)
print(SOC_10_20cm)
sink()
#-------------------------------------------------------------------------------
# Create Shoot Biomass Figures
soc_gg = ggplot(soilC_dt[Depth == '0to10'],aes(x = Trt, y = soc.Mg.ha, fill = Trt, group = Trt)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot(alpha = 0.9) +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  ylim(0.0,35) +
  facet_grid(~ Year, labeller = label_parsed) +
  scale_fill_viridis(discrete = TRUE, option = "D", begin = 0.6, end = 0.2) +
  ylab(expression("Soil Organic Carbon (Mg C ha"^-1*")")) +
  xlab("Treatment") +
  theme_bw() +
  theme(text=element_text(size=7, color = 'black'),
        axis.title = element_text(size = 8, color = 'black'),
        strip.text.x = element_text(size = 7, color = 'black'),
        axis.text=element_text(size=7, color = 'black'),
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        # axis.title.y =  element_blank(),
        # axis.title.x = element_blank(),
        legend.title=element_text(size = 7, color = 'black'),
        legend.position = "none",
        plot.title.position = "plot",  # This moves the title to align with plot edge
        plot.title = element_text(
          hjust = -0.005,  # Slight adjustment left of the plot
          vjust = -0.5,   # Slight adjustment above the plot
          size = 7       # Match your other text size if needed
        )) +
  guides(fill=guide_legend(title="Treatment"))
soc_gg

# add P-values
text = data.table(Trt = 'Compost', soc.Mg.ha = 32, lab = "P = 0.001",
                  Year  = factor(c("2018"),levels = c("2018","2020")))
soc_gg = soc_gg + geom_text(data = text, label = "P = 0.007", 
                                        nudge_x = 0.5, size = 3)
text = data.table(Trt = 'Compost', soc.Mg.ha = 32, lab = "P = 0.001",
                  Year  = factor(c("2020"),levels = c("2018","2020")))
soc_gg = soc_gg + geom_text(data = text, label = "P < 0.001", 
                                        nudge_x = 0.5, size = 3)
soc_gg

ggsave("FigureS2.tiff", plot = soc_gg, path = figures_path, width = 88, height = 88, units = "mm", dpi = 600)
