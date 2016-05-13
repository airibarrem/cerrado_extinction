#rm(list=ls())
require(raster)
require(snowfall)
source('./cerrado_fun.R')

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


# Vegetation on 2010
veg.shp = shapefile('./shapes/2_Cerrado_in_2012.shp')

# Areas to be restored by 2050
rest.shp = shapefile('./shapes/7_restoration.shp')

# Reprojecting the vegetation shapes
veg.shp = spTransform(veg.shp,CRS(proj4string(bgd)))
rest.shp = spTransform(rest.shp,CRS(proj4string(bgd)))

veg.ras = rasterize(veg.shp[1],bgd)
rest.ras = rasterize(rest.shp[1],bgd)

veg.ras = veg.ras/veg.ras
veg.ras[is.na(veg.ras)]=0
veg.ras = veg.ras+bgd

rest.ras = rest.ras/rest.ras
rest.ras[is.na(rest.ras)]=0
rest.ras=rest.ras+bgd

# Computing vegetation raster by 2050 under restoration scenario
veg.rest = veg.ras + rest.ras
veg.rest = veg.rest / veg.rest


# Land-use raster
OA.now = rcompose('./LU/uso_da_terra_culturas2safra.tif',bgd)

# native vegetation mask (according to classes in original land-use rasters)
OA.nat.now = (OA.now %in% c(5,6,7,8)) + bgd

rest.run = hab.calc(flora.files, OA.nat.now, veg.rest, bgd, binary.px=T,
                   sp.dir=flora.dir, string.rem='.asc', CSV.name='Flora_OA_Rest',
                   area.only=F, print.CSV=T, sf.on=T, cores=7)


