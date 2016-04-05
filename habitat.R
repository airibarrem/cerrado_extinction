#rm(list=ls())
require(raster)
require(snowfall)
source('./cerrado_fun.R')

start = Sys.time()

###


# Getting the list of names of all sp.dist files to be included
flora.dir = './Sp_maps/Flora_Cerrado/'
fauna.dir = './Sp_maps/Fauna_Cerrado/'
birds.dir = './Sp_maps/Aves_Cerrado/'
mmals.dir = './Sp_maps/Mammals_Cerrado/MML/'

flora.files = dir(flora.dir)
fauna.files = dir(fauna.dir)
birds.files = dir(birds.dir)
mmals.files = dir(mmals.dir)

# Keeping only .asc filenames
flora.files = flora.files[grepl('.asc',flora.files)]
fauna.files = fauna.files[grepl('.asc',fauna.files)]
birds.files = birds.files[grepl('.asc',birds.files)]
mmals.files = mmals.files[grepl('.asc',mmals.files)]

# Removing future distributions from mammals filenames:
mmals.files = mmals.files[!grepl('_mml_f.asc',mmals.files)]


###


# Rasterizing the background shapefile Using first file of sp.dist. as template
bgd = rasterize(shapefile('./biome/Cerrado.shp'),
                raster(paste0(mmals.dir,mmals.files[1])))

# preparing bgd to be used as addition mask (add to a raster to apply the mask)
bgd = bgd/bgd
bgd[bgd>0]=0


###


# REDD-PAC land-use shapefile
BB.shp = shapefile('./LU_REDD_PAC/REED_PAC_BAU.shp')

# Reprojecting LU shapefiles to match background of sp.dist rasters
BB.shp = spTransform(BB.shp,CRS(proj4string(bgd)))

# Names of classes (columns) in LU.shp
BB.cls = names(BB.shp)

# Native vegetation shapefiles (now = 2010, bau = 2050)
BB.nat.shp.now = BB.shp[BB.cls %in% c('NatLnd2010','PRIFOR2010','ForReg2010')]
BB.nat.shp.bau = BB.shp[BB.cls %in% c('NatLnd2050','PRIFOR2050','ForReg2050')]

# Converting native vegetation shapefiles to a combined raster of veg. cover
BB.nat.shp.now = rasterize(BB.nat.shp.now,bgd)
BB.nat.shp.bau = rasterize(BB.nat.shp.bau,bgd)

# Each subs command puts in each pixel its corresponding value computed by rasterize
# Sum over each class in nat.shp corresponding to column in @data@attributes
BB.nat.now = subs(BB.nat.shp.now,BB.nat.shp.now@data@attributes[[1]][,c(1,2)]) +
  subs(BB.nat.shp.now,BB.nat.shp.now@data@attributes[[1]][,c(1,3)]) +
  subs(BB.nat.shp.now,BB.nat.shp.now@data@attributes[[1]][,c(1,4)]) + bgd
BB.nat.bau = subs(BB.nat.shp.bau,BB.nat.shp.bau@data@attributes[[1]][,c(1,2)]) +
  subs(BB.nat.shp.bau,BB.nat.shp.bau@data@attributes[[1]][,c(1,3)]) +
  subs(BB.nat.shp.bau,BB.nat.shp.bau@data@attributes[[1]][,c(1,4)]) + bgd

# Converting area_size (in ha) shapefile to a raster, with area data in all px
size.map = rasterize(BB.shp[30],bgd)
size.map = subs(size.map, size.map@data@attributes[[1]][,1:2]) + bgd

# Converting pixel size in native vegetation rasters
# Multiplicative factor 1000 due to unit differences in original shape file
# now each pixel has a value of veg. area (ha) corresponding to its own area
BB.nat.now = (BB.nat.now * 1000) * ((area(bgd) * 100) / size.map)
BB.nat.bau = (BB.nat.bau * 1000) * ((area(bgd) * 100) / size.map)


###


# Land-use rasters
OA.now = rcompose('./LU/uso_da_terra_culturas2safra.tif',bgd)
OA.bau = rcompose('./LU/Brasil2050.tif',bgd)

# native vegetation mask (according to classes in original land-use rasters)
OA.nat.now = (OA.now %in% c(5,6,7,8)) + bgd
OA.nat.bau = (OA.bau %in% c(5,6,7,8)) + bgd


###


# Preparing input lists

dir.list = list('Flora'=flora.dir, 'Fauna'=fauna.dir, 'Birds'=birds.dir,
                'Mammals'=mmals.dir)

files.list = list('Flora'=flora.files, 'Fauna'=fauna.files, 'Birds'=birds.files,
                  'Mammals'=mmals.files)

bin.opt.list = list('BRABIOM'=F, 'Otimizagro'=T)

str.opt.list = list('Flora'='.asc', 'Fauna'='.asc', 'Birds'='_p.asc',
                    'Mammals'='_mml_p.asc')

nat.now.list = list('BRABIOM'=BB.nat.now, 'Otimizagro'=OA.nat.now)
nat.bau.list = list('BRABIOM'=BB.nat.bau, 'Otimizagro'=OA.nat.bau)


###


# Initiating final list of data.frames
df.list = list('BRABIOM'=0, 'Otimizagro'=0)

# Combining all groups vs LU maps
for (gr in 1:length(dir.list)){
  for (lu in 1:length(nat.now.list)){
    
    # Printing ellapsed time
    print(Sys.time()-start)
    
    # Setting inputs for the current run
    sp.files = files.list[[gr]]
    sp.dir = dir.list[[gr]]
    nat.now = nat.now.list[[lu]]
    nat.bau = nat.bau.list[[lu]]
    bin.opt = bin.opt.list[[lu]]
    str.opt = str.opt.list[[gr]]
    filename = paste0(names(files.list)[[gr]],'_',names(nat.bau.list)[[lu]])
    print(filename)
    
    # Computing habitats for the current run
    hab.run = hab.calc(sp.files, nat.now, nat.bau, binary.px=bin.opt,
                       sp.dir=sp.dir, string.rem=str.opt, CSV.name = filename,
                       area.only=T, print.CSV=T, sf.on=T, cores=11)
    
    # Registering group of the current run
    hab.run$group = names(files.list)[[gr]]
    
    # Combining habitat results from current run to final data.frame
    if (length(df.list[[lu]]==1)) {df.list[[lu]]=hab.run} else {df.list[[lu]]=rbind(df.list[[lu]],hab.run)}
  }
}

# Organinzing final data.frame by species name
hab.df = hab.df[sort(hab.df$Sp),]

# Outputing final data.frame
write.csv(hab.df, file=paste0('Combined_results_',Sys.Date()) )