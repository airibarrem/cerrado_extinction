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


# Land-use rasters
LU.now = rcompose('./LU/uso_da_terra_culturas2safra.tif',bgd)
LU.bau = rcompose('./LU/Brasil2050.tif',bgd)

# native vegetation mask (according to classes in original land-use rasters)
nat.now = (LU.now %in% c(6,7,8)) + bgd
nat.bau = (LU.bau %in% c(6,7,8)) + bgd


# Applying intersect.raster and cond.area to all sp.dist files
sfInit(parallel = T, cpus = 3, type = 'SOCK')
  sfExportAll()
  sfLibrary(raster)
  # sp. potential area maps today
  pot.maps = sfLapply(paste0('./Flora_Cerrado/',sp.files),rcompose,bgd)
  # sp. habitats today
  hab.maps = sfLapply(paste0('./Flora_Cerrado/',sp.files),intersect.raster,nat.now,bgd)
  #sp. habitats under future LUC
  bau.maps = sfLapply(paste0('./Flora_Cerrado/',sp.files),intersect.raster,nat.bau,bgd)
  # List of all the maps for today
  maps.list = c(pot.maps,hab.maps,bau.maps)
  # List of the areas (in ha) corresponding to the list of maps above
  pot.areas = sfLapply(pot.maps,function(x){cond.area(nat.now,x)})
  hab.areas = sfLapply(hab.maps,function(x){cond.area(nat.now,x)})
  bau.areas = sfLapply(bau.maps,function(x){cond.area(nat.bau,x)})
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
write.csv(now.df, file='Flora_otimizagro')
write.csv2(now.df, file='Flora_otimizagro_2')

write.csv(teste.df, file='Flora_otimizagro_perc')
write.csv2(teste.df, file='Flora_otimizagro_perc_2')
