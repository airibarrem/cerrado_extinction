#rm(list=ls())
require(raster)
require(snowfall)

heaviside = function(x){stepfun(0,c(0,1))(x)}

rcompose = function(x, bgx){
  # read input raster, if x is filename
  if (is.character(x)){x = raster(x)}
  # projects input raster into background raster
  if (is.na(crs(x))) {crs(x)=crs(bgx)}
  if (grepl(crs(x),crs(bgx))) {x=resample(x,bgx,method='ngb')} else {x=projectRaster(x,bgx,method='ngb')}
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

# Computes habitat rasters and area (in ha)
hab.calc = function(sp.filenames, nat.now, nat.bau, bgd, binary.px=T, area.only=F,
                    sp.dir=getwd(), string.rem='.asc', CSV.name = 'hab_areas',
                    print.CSV=T, sf.on=T, cores=11){
  # INPUTS:
  # sp.filenames    list of filenames of the sp. dist. rasters of be included
  # nat.now         raster corresponding to natural vegetation in the present
  # nat.bau         raster corresponding to natural vegetation in the BAU
  # sp.dir          character line containing filepath to sp.files
  # string.rem      string to be removed from sp. filenames in final print
  # binary.px       px in nat. rasters: 0/1 presence (T), 0-1 percent cover (F)
  # area.only       If T, returns only the aggregated area data.frame
  # print.CSV       logical. prints csv with areas in ha when T
  # CSV.name        filename for the CSV to be printed
  # sf.on           logical. uses snowfall parallel computing when T
  # cores           number of cores to be used in snowfall calls
  #
  # Applying intersect.raster and cond.area to all sp.dist files
  sfInit(parallel = sf.on, cpus = cores, type = 'SOCK')
  sfExport('sp.dir','sp.files','rcompose','bgd','intersect.raster','nat.now',
           'nat.bau','cond.area')
  sfLibrary(raster)
  # sp. potential area maps today
  pot.maps = sfLapply(paste0(sp.dir,sp.files),rcompose,bgd)
  # sp. habitats today
  hab.maps = sfLapply(paste0(sp.dir,sp.files),intersect.raster,nat.now,bgd)
  #sp. habitats under future LUC
  bau.maps = sfLapply(paste0(sp.dir,sp.files),intersect.raster,nat.bau,bgd)
  # List of all the maps for today
  maps.list = c(pot.maps,hab.maps,bau.maps)
  # List of the areas (in ha) corresponding to the list of maps above
  pot.areas = sfLapply(pot.maps,function(x){cond.area(nat.now,x)})
  if (binary.px){
    hab.areas = sfLapply(hab.maps,function(x){cond.area(nat.now,x)})
    bau.areas = sfLapply(bau.maps,function(x){cond.area(nat.bau,x)})
  } else {
    hab.areas = sfLapply(hab.maps,function(x){sum(values(x),na.rm=T)})
    bau.areas = sfLapply(bau.maps,function(x){sum(values(x),na.rm=T)})
  }
  sfStop()
  
  # Computing area data.frames
  res.df = data.frame(Sp=sub("_",' ',sub(string.rem,'',sp.files)))
  res.df$Pot_Area = as.numeric(pot.areas)
  res.df$Hab_Area = as.numeric(hab.areas)
  res.df$BAU_Area = as.numeric(bau.areas)
  
  res.df$Hab_perc = res.df$Hab_Area / res.df$Pot_Area
  res.df$BAU_perc = res.df$BAU_Area / res.df$Pot_Area
  
  res.df$Ext_015 = 1 - (res.df$BAU_perc^0.15)
  res.df$Ext_025 = 1 - (res.df$BAU_perc^0.25)
  res.df$Ext_035 = 1 - (res.df$BAU_perc^0.35)
  
  res.df = res.df[order(res.df$Sp),]
  
  # Outputting area data.frames
  write.csv(res.df, file=CSV.name)
  
  # Packing results
  if (area.only) {res=res.df} else {res=list('data'=res.df, 'maps'=maps.list)}
  
  return(res)
}
