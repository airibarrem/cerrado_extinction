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
sp.files = dir('./Mammals_Cerrado/MML')

# Splitting files in (p)resent and (f)uture maps
sp.now = sp.files[grep('_p.asc',sp.files)]
sp.bau = sp.files[grep('_f.asc',sp.files)]


# Rasterizing the background shapefile Using first file of sp.dist. as template
bgd = rasterize(shapefile('./biome/Cerrado.shp'),
                raster(paste0('./Mammals_Cerrado/MML/',sp.files[1])))

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
  pot.now = sfLapply(paste0('./Mammals_Cerrado/MML/',sp.now),rcompose,bgd)
  # sp. habitats today
  hab.now = sfLapply(paste0('./Mammals_Cerrado/MML/',sp.now),intersect.raster,nat.now,bgd)
  #sp. habitats under future LUC
  luc.now = sfLapply(paste0('./Mammals_Cerrado/MML/',sp.now),intersect.raster,nat.bau,bgd)
  # List of all the maps for today
  maps.now = c(pot.now,hab.now)
  # List of the areas (in ha) corresponding to the list of maps above
  areas.now = sfLapply(maps.now,function(x){cond.area(nat.now,x)})
  areas.now = c(areas.now,sfLapply(luc.now,function(x){cond.area(nat.bau,x)}))
  #
  # Same procedure, but for projected maps (BAU 2050)
  pot.bau = sfLapply(paste0('./Mammals_Cerrado/MML/',sp.bau),rcompose,bgd)
  hab.bau = sfLapply(paste0('./Mammals_Cerrado/MML/',sp.bau),intersect.raster,nat.bau,bgd)
  # Effect of climate change under a zero deforestation (LNAE) scenario
  cc.lnae = sfLapply(paste0('./Mammals_Cerrado/MML/',sp.bau),intersect.raster,nat.now,bgd)
  maps.bau = c(pot.bau,hab.bau)
  areas.bau = sfLapply(maps.bau,function(x){cond.area(nat.bau,x)})
  areas.bau = c(areas.bau,sfLapply(cc.lnae,function(x){cond.area(nat.bau,x)}))
sfStop()

# Computing area data.frames
now.df = data.frame(Sp=sub("_",' ',sub('_mml_p.asc','',sp.now)))
now.df$Pot_Area = unlist(areas.now)[1:length(sp.now)]
now.df$Hab_Area = unlist(areas.now)[(1+length(sp.now)):(2*length(sp.now))]
now.df$LUC_Area = unlist(areas.now)[(1+2*length(sp.now)):(3*length(sp.now))]
now.df = now.df[order(now.df$Sp),]

bau.df = data.frame(Sp=sub("_",' ',sub('_mml_f.asc','',sp.bau)))
bau.df$Pot_Area = unlist(areas.bau)[1:length(sp.bau)]
bau.df$Hab_Area = unlist(areas.bau)[(1+length(sp.bau)):(2*length(sp.bau))]
bau.df$LNAE_Area = unlist(areas.bau)[(1+2*length(sp.bau)):(3*length(sp.bau))]
bau.df = bau.df[order(bau.df$Sp),]

comb.df = data.frame(Sp=now.df[now.df$Sp %in% bau.df$Sp,1])
comb.df$Now_Pot = now.df[now.df$Sp %in% bau.df$Sp,2]
comb.df$Now_Hab = now.df[now.df$Sp %in% bau.df$Sp,3]
comb.df$Int_LUC = now.df[now.df$Sp %in% bau.df$Sp,4]
comb.df$CC_LNAE = bau.df[,4]
comb.df$BAU_Pot = bau.df[,2]
comb.df$BAU_Hab = bau.df[,3]

teste.df = comb.df
teste.df[,3:6] = comb.df[,3:6] / comb.df[,2]
teste.df[,7] =  comb.df[,7] / comb.df[,6]

# Outputting area data.frames
write.csv(now.df, file='Mammals_Present')
write.csv(now.df, file='Mammals_BAU')
write.csv(comb.df, file='Mammals_combined')
write.csv(teste.df, file='Mammals_combined_perc')

write.csv2(now.df, file='Mammals_Present_2')
write.csv2(now.df, file='Mammals_BAU_2')
write.csv2(comb.df, file='Mammals_combined_2')
write.csv2(teste.df, file='Mammals_combined_perc_2')
