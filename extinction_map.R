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
OA.nat.now = (OA.now %in% c(5,6,7,8)) + bgd
OA.nat.bau = (OA.bau %in% c(5,6,7,8)) + bgd

writeRaster(OA.nat.now,filename='Otimizagro_natveg_2010.tif',overwrite=T)
writeRaster(OA.nat.bau,filename='Otimizagro_natveg_2050.tif',overwrite=T)

jpeg("./Otimizagro_natveg.jpg", res=600, width=15, height=9, unit='cm', pointsize=8)
spplot(stack(OA.nat.now,OA.nat.bau,OA.nat.now-OA.nat.bau), names.attr=c(2010,2050,'Deforestation'),
       col.regions=rainbow(18,start=0.1, end=0.3))
dev.off()

hab.run = hab.calc(flora.files, OA.nat.now, OA.nat.bau, binary.px=T,
                   sp.dir=flora.dir, string.rem='.asc', CSV.name='Flora_Otimizagro',
                   area.only=F, print.CSV=T, sf.on=T, cores=7)

pot.maps = unlist(hab.run[[2]][1:length(ext.data[,5])])
hab.maps = unlist(hab.run[[2]][(length(ext.data[,5])+1):(2*length(ext.data[,5]))])
bau.maps = unlist(hab.run[[2]][(2*length(ext.data[,5])+1):(3*length(ext.data[,5]))])

ext.pot.maps = pot.maps
ext.hab.maps = hab.maps
ext.bau.maps = bau.maps

for (i in 1:length(ext.data[,5])) {ext.pot.maps[[i]] = ext.data[i,5] * pot.maps[[i]]}
for (i in 1:length(ext.data[,5])) {ext.hab.maps[[i]] = ext.data[i,5] * hab.maps[[i]]}
for (i in 1:length(ext.data[,5])) {ext.bau.maps[[i]] = ext.data[i,5] * bau.maps[[i]]}

rich.pot = Reduce('+',pot.maps)
rich.hab = Reduce('+',hab.maps)
rich.bau = Reduce('+',bau.maps)

writeRaster(rich.pot, filename='richness.pot.tif',overwrite=T)
writeRaster(rich.hab, filename='richness.hab.tif',overwrite=T)
writeRaster(rich.bau, filename='richness.bau.tif',overwrite=T)


final.pot.map = Reduce('+',ext.pot.maps)
final.hab.map = Reduce('+',ext.hab.maps)
final.bau.map = Reduce('+',ext.bau.maps)

writeRaster(final.pot.map, filename='final.pot.map.tif',overwrite=T)
writeRaster(final.hab.map, filename='final.hab.map.tif',overwrite=T)
writeRaster(final.bau.map, filename='final.bau.map.tif',overwrite=T)

jpeg("./extinction_map_pot.jpg", res=600, width=9, height=9, unit='cm', pointsize=8)
spplot(final.pot.map, names.attr='Aggregated extinction probabilities 2050',
       col.regions=colorRampPalette(colors=c('ivory2','red4','red'),interpolate='spline'))
dev.off()

jpeg("./extinction_map_hab.jpg", res=600, width=9, height=9, unit='cm', pointsize=8)
spplot(final.hab.map, names.attr='Aggregated extinction probabilities 2050',
       col.regions=colorRampPalette(colors=c('ivory2','red4','red'),interpolate='spline'))
dev.off()

jpeg("./extinction_map_bau.jpg", res=600, width=9, height=9, unit='cm', pointsize=8)
spplot(final.bau.map, names.attr='Aggregated extinction probabilities 2050',
       col.regions=colorRampPalette(colors=c('ivory2','red4','red'),interpolate='spline'))
dev.off()

jpeg("./extinction_map_comb.jpg", res=600, width=15, height=9, unit='cm', pointsize=8)
spplot(stack(final.pot.map,final.hab.map,final.bau.map), names.attr=c('Potential','2010','2050'),
       col.regions=colorRampPalette(colors=c('ivory2','red4','red'),interpolate='spline'))
dev.off()

jpeg("./richness_map_comb.jpg", res=600, width=15, height=9, unit='cm', pointsize=8)
spplot(stack(rich.pot,rich.hab,rich.bau), names.attr=c('Potential','2010','2050'),
       col.regions=colorRampPalette(colors=c('ivory2','green4','green'),interpolate='spline'))
dev.off()
