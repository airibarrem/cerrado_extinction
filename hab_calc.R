# Computes habitat rasters and area (in ha)
hab.calc = function(sp.files, LU.file, bgd.file, nat.classes = c(1),
                    sp.dir=getwd(), LU.dir=getwd(), bgd.dir=getwd(),
                    string.rem='', print.CSV=T, sf.on=T, cores=3){
  # INPUTS:
  # sp.files    list of filenames of the sp. dist. rasters of be included
  # LU.file     filename with the map (.shp or .asc) of the land-use to be used
  # bgd.file    filename with the map (.shp or .asc) of contour of the region
  # nat.classes list with the class IDs in LU.file corresponding to natural veg.
  # sp.dir      character line containing filepath to sp.files
  # LU.dir      character line containing filepath to LI.file
  # bgd.dir     character line containing filepath to bgd.files
  # string.rem  string to be removed from sp. filenames in final print
  # print.CSV   logical. prints csv with areas in ha when T
  # sf.on       logical. uses snowfall parallel computing when T
  # cores       number of cores to be used in snowfall calls
  #
  #
  library('raster','snowfall')
  #
  # Preparing background raster
  bgd = ifelse(grepl('.shp',bgd.file),
               rasterize(shapefile(paste0(bgd.dir,bgd.file)),
                         raster(paste0(sp.dir,sp.files[1]))),
               rcompose(paste0(bgd.dir,bgd.file)))
  # preparing bgd to be used as addition mask (add to a raster to apply mask)
  bgd = bgd/bgd
  bgd[bgd>0]=0
  # Land-use raster
  LU.ras = ifelse(grepl('.shp',LU.file),
                  rcompose(rasterize(shapefile(paste0(LU.dir,LU.file)),
                                     raster(paste0(sp.dir,sp.files[1]))),bgd),
                  rcompose(paste0(LU.dir,LU.file),bgd))
  # native vegetation mask (according to classes in original land-use rasters)
  nat.ras = (LU.ras %in% nat.classes) + bgd
  # Starts snowfall cluster, and exports packages and variables to all cores
  sfInit(parallel = sf.on, cpus = cores, type = 'SOCK')
  sfExportAll()
  sfLibrary(raster)
  # list with rasterLayers for all species in sp.files
  sp.cube = sfLapply(paste0(sp.dir,sp.files),rcompose,bgd)
  # list with rasterLayers with sp. habitats (combining sp.cube w nat.veg map)
  hab.cube = sfLapply(sp.cube, intersect.area,nat.ras,bgd)
  #sp. habitats under future LUC
  # List of all the maps
  maps.cube = c(sp.cube,hab.cube)
  # List of the areas (in ha) corresponding to the list of maps above
  areas.now = sfLapply(maps.cube,function(x){cond.area(nat.ras,x)})
  sfStop()
  # Organizing output data.frame
  res.df = data.frame(Sp=sub("_",' ',sub(string.rem,'',sp.now)))
  res.df$Pot_Area = unlist(areas.now)[1:length(sp.now)]
  res.df$Hab_Area = unlist(areas.now)[(1+length(sp.now)):(2*length(sp.now))]
  res.df = now.df[order(now.df$Sp),]
  # Printing output CSV
  if(print.CSV){write.csv(res.df, file=paste0('hab_calc_results_',Sys.date()) )}
  return(list('Areas_table'=res.df,'Sp._rasters'=sp.cube, 'Habitat_rasters'=hab.cube))
}

hab.calc(sp.bau,'Brasil2050.tif','Cerrado.shp',nat.classes=c(6,7,8),
         sp.dir='./Mammals_Cerrado/MML',LU.dir='./LU/',bgd.dir='.biome',
         string.rem='_mml_f.asc')
