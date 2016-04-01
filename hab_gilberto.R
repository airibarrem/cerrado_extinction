rm(list=ls())
require(raster)
require(snowfall)

heaviside = function(x){stepfun(0,c(0,1))(x)}

rcompose = function(x, bgx){
  # read input raster, if x is filename
  if (is.character(x)){x = raster(x)}
  # projects input raster into background raster
  if(!is.na(crs(x))){x = projectRaster(x,bgx,method='ngb')}
  # remove NA values from input raster
  x[is.na(x)]=0
  # mask input raster using the background one
  x = x + bgx
  return(x)
}

cond.area = function(x, xcond, convert.ha=TRUE, px.val=F){
  # Returns an approximate area of pixels in "x" that satisfy "xcond" condition.
  # "xcond" must be a logical raster (px values either 0 or 1).
  zone.count = 0
  zone.res = zonal(x,xcond,fun=function(x, na.rm){ if(na.rm){length(na.omit(x))}
                                                   else{length(x)}})
  if(length(zone.res[,1]) > 1){zone.count = zone.res[[2,2]]}
  px.a = mean(values(area(x)))
  res = px.a*zone.count
  if (convert.ha==TRUE){res = res*100}
  if (px.val == T){
    print(paste("Px count:",zone.count))
    res = zone.count}
  return(res)
}

# Intersects x with y, and masks the result with bgx
intersect.raster = function(x, y, bgx){
  if(is.character(x)){x=rcompose(x,bgx)}
  if(is.character(y)){y=rcompose(y,bgx)}
  return(rcompose(x*y,bgx))}


###


# Getting the list of names of all sp.dist files to be included
sp.files = dir('./Flora_Cerrado/')

# Rasterizing the background shapefile Using first file of sp.dist. as template
bgd = rasterize(shapefile('./biome/Cerrado.shp'),
                raster(paste0('./Flora_Cerrado/',sp.files[1])))

# preparing bgd to be used as addition mask (add to a raster to apply the mask)
bgd = bgd/bgd
bgd[bgd>0]=0


# REDD-PAC land-use shapefile
LU.shp = shapefile('./LU_REDD_PAC/REED_PAC_BAU.shp')

# Reprojecting LU shapefiles to match background of sp.dist rasters
LU.shp = spTransform(LU.shp,CRS(proj4string(bgd)))

# Names of classes (columns) in LU.shp
LU.classes = names(LU.shp)


# Native vegetation shapefiles (now = 2010, bau = 2050)
nat.shp.now = LU.shp[LU.classes %in% c('NatLnd2010','PRIFOR2010','ForReg2010')]

nat.shp.bau = LU.shp[LU.classes %in% c('NatLnd2050','PRIFOR2050','ForReg2050')]

# Converting native vegetation shapefiles to a combined raster of veg. cover
nat.shp.now = rasterize(nat.shp.now,bgd)

# Each subs command puts in each pixel its corresponding value computed by rasterize
# Sum over each class in nat.shp corresponding to column in @data@attributes
nat.now = subs(nat.shp.now,nat.shp.now@data@attributes[[1]][,c(1,2)]) +
  subs(nat.shp.now,nat.shp.now@data@attributes[[1]][,c(1,3)]) +
  subs(nat.shp.now,nat.shp.now@data@attributes[[1]][,c(1,4)]) + bgd

# Same procedure for BAU-projection classes
nat.shp.bau = rasterize(nat.shp.bau,bgd)

nat.bau = subs(nat.shp.bau,nat.shp.bau@data@attributes[[1]][,c(1,2)]) +
  subs(nat.shp.bau,nat.shp.bau@data@attributes[[1]][,c(1,3)]) +
  subs(nat.shp.bau,nat.shp.bau@data@attributes[[1]][,c(1,4)]) + bgd

# Converting area_size (in ha) shapefile to a raster, with data in all pixels
size.map = rasterize(LU.shp[30],bgd)
size.map = subs(size.map, size.map@data@attributes[[1]][,1:2]) + bgd

# Converting pixel size in native vegetation rasters
# Multiplicative factor 1000 due to unit differences in original shape file
# now each pixel has a value of veg. area (ha) corresponding to its own area
nat.now = (nat.now * 1000) * ((area(bgd) * 100) / size.map)
nat.bau = (nat.bau * 1000) * ((area(bgd) * 100) / size.map)


# Applying intersect.raster and cond.area to all sp.dist files
sfInit(parallel = T, cpus = 7, type = 'SOCK')
sfExportAll()
sfLibrary(raster)
# sp. potential area maps today (each pixel contains its area in ha)
pot.maps = sfLapply(paste0('./Flora_Cerrado/',sp.files),
                    function(x){rcompose(x,bgd)*area(bgd)*100})
# sp. habitats today
hab.maps = sfLapply(paste0('./Flora_Cerrado/',sp.files),
                    intersect.raster,nat.now,bgd)
#sp. habitats under future LUC
bau.maps = sfLapply(paste0('./Flora_Cerrado/',sp.files),
                    intersect.raster,nat.bau,bgd)
# List of all the maps for today
maps.list = c(pot.maps,hab.maps,bau.maps)
# List of the areas (in ha) corresponding to the list of maps above
pot.areas = sfLapply(pot.maps,function(x){sum(values(x),na.rm=T)})
hab.areas = sfLapply(hab.maps,function(x){sum(values(x),na.rm=T)})
bau.areas = sfLapply(bau.maps,function(x){sum(values(x),na.rm=T)})
sfStop()

# Computing area data.frames
now.df = data.frame(Sp=sub("_",' ',sub('.asc','',sp.files)))
now.df$Pot_Area = as.numeric(pot.areas)
now.df$Hab_Area = as.numeric(hab.areas)
now.df$BAU_Area = as.numeric(bau.areas)
now.df = now.df[order(now.df$Sp),]

teste.df = now.df
teste.df[,3:4] = now.df[,3:4] / as.numeric(now.df[,2])


# Outputting area data.frames
write.csv(now.df, file='Flora_BRABIOM')
write.csv2(now.df, file='Flora_BRABIOM_2')

write.csv(teste.df, file='Flora_BRABIOM_perc')
write.csv2(teste.df, file='Flora_BRABIOM_perc_2')

# Exporting native vegetation maps
jpeg("./BRABIOM_natveg.jpg", res=600, width=9, height=9, unit='cm', pointsize=8)
spplot(stack(nat.now,nat.bau), names.attr=c(2010,2050),
       col.regions=rainbow(18,start=0.1, end=0.3))
dev.off()
