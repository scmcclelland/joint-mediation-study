# Map figure 
# Author: Shelby C McClelland
# Created:     23 February 2025
# Last Update: 23 February 2025
# Description: This file creates components of map for Figure S1.
#-------------------------------------------------------------------------------
library(cowplot)
library(data.table)
library(ggplot2)
library(rstudioapi)
library(usmap)
base_path    = dirname(getActiveDocumentContext()$path)
figures_path = paste(base_path, 'figures', sep = '/')
#-------------------------------------------------------------------------------
# us map
#-------------------------------------------------------------------------------
dt = usmapdata::fips_data(regions = "counties")
dt = setDT(as.data.table(dt))

co_dat      = data.table(state = c("CO"), f = c(1)) # create fill data
st_dat      = data.table(state = unique(dt[!abbr == "CO",abbr]), f = 0)
#rbind
st_dat      = rbind(st_dat, co_dat)
setorder(st_dat, state)
base_map = plot_usmap(data = st_dat, 
                      # fill = "white",
                      values = "f",
                      color = "gray60",
                      regions = "states",
           exclude = c('Alaska', 'Hawaii')) +
    scale_fill_continuous(
      low = "white", high = "gray40", name = "f"
    ) # add fill color
base_map = base_map + theme_bw() + theme(legend.position = 'none',
                                         axis.text = element_text(size = 10))
base_map
#-------------------------------------------------------------------------------
# co map
#-------------------------------------------------------------------------------
dt     = dt[full == 'Colorado'] #fips = 08069 (Larimer County)
lar_dat = data.table(county = c('Larimer County'), fips = 08069, f = c(1))
cou_dat = data.table(county = unique(dt[!county == 'Larimer County', county]),
                     fips = unique(dt[!county == 'Larimer County', fips]),
                    f = c(0))
#rbind
cou_dat = rbind(cou_dat, lar_dat)
setorder(cou_dat, county)

# fill in Larimer county
co_map = plot_usmap(data = cou_dat, 
                      # fill = "white",
                      values  = "f",
                      color   = "gray60",
                      regions = "counties",
                      include = "08") +
  scale_fill_continuous(
    low = "white", high = "gray40", name = "f"
  ) # add fill color
co_map = co_map + theme_bw() + theme(legend.position = 'none',
                                     axis.text = element_text(size = 10))
co_map

combined_maps = ggdraw() +
  draw_plot(base_map) +
  draw_plot(co_map,
            x = 0.60,      # x position of inset
            y = 0.11,      # y position of inset
            width = 0.4,   # height of inset
            height = 0.4)
combined_maps
ggsave("FigureS1.tiff", plot = combined_maps, path = figures_path, width = 180, height = 90, units = "mm", dpi = 600)
