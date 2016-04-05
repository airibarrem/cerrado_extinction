#rm(list=ls())
require(raster)
require(snowfall)
source('./cerrado_fun.R')

ext.data = read.csv('./Flora_soaresfilhox.csv')[,2:6]

# Getting the list of names of all sp.dist files to be included
flora.dir = './Sp_maps/Flora_Cerrado/'

flora.files = dir(flora.dir)

# Keeping only .asc filenames
flora.files = flora.files[grepl('.asc',flora.files)]

# Rasterizing the background shapefile Using first file of sp.dist. as template
bgd = rasterize(shapefile('./biome/Cerrado.shp'),
                raster(paste0(flora.dir,flora.files[1])))

# preparing bgd to be used as addition mask (add to a raster to apply the mask)
bgd = bgd/bgd
bgd[bgd>0]=0

# Land-use rasters
OA.now = rcompose('./LU/uso_da_terra_culturas2safra.tif',bgd)
OA.bau = rcompose('./LU/Brasil2050.tif',bgd)

# native vegetation mask (according to classes in original land-use rasters)
OA.nat.now = (OA.now %in% c(6,7,8)) + bgd
OA.nat.bau = (OA.bau %in% c(6,7,8)) + bgd

hab.run = hab.calc(flora.files, OA.nat.now, OA.nat.bau, binary.px=T,
                   sp.dir=flora.dir, string.rem='.asc',
                   area.only=F, print.CSV=F, sf.on=T, cores=7)

bau.maps = unlist(hab.run[[2]][(2*length(ext.data[,5]))+1:(3*length(ext.data[,5]))])

ext.maps = bau.maps

for (i in 1:length(ext.data[,5])) {ext.maps[[i]] = ext.data[i,5] * bau.maps[[i]]}

final.map = Reduce('+',ext.maps)

jpeg("./extinction_map.jpg", res=600, width=9, height=9, unit='cm', pointsize=8)
spplot(final.map, names.attr='Aggregated extinction probabilities 2050',
       col.regions=colorRampPalette(colors=c('seashell3','red4','red')),add=T,alpha=0.1)
dev.off()
