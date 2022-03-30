rm(list = ls())

require(sf)
require(terra)


countries <- st_read("../../Datasets/Country_shapefiles/all_countries.shp")


sao_tome <- countries %>% dplyr::filter(grepl(iso_short, pattern = "sao tome", ignore.case = TRUE))
sao_tome <- as_Spatial(sao_tome$geometry)

sao_tome_extent <- extent(6.4,6.9,0.048,1)
principe_extent <- extent(7,7.45,1,1.699121)

principe <- crop(sao_tome,principe_extent)
crs(principe) <- "+proj=longlat +datum=WGS84 +no_defs"

sao_tome <- crop(sao_tome,sao_tome_extent)
crs(sao_tome) <- "+proj=longlat +datum=WGS84 +no_defs"


comoros <- countries %>% dplyr::filter(grepl(iso_short, pattern = "comoros", ignore.case = TRUE))
comoros <- as_Spatial(comoros$geometry)
comoros_extent <- extent(43.22,43.7,-12,-11.36)
comoros <- crop(comoros,comoros_extent)
crs(comoros) <- "+proj=longlat +datum=WGS84 +no_defs"

puerto_rico <- countries %>% dplyr::filter(grepl(iso_short, pattern = "puerto rico", ignore.case = TRUE))
puerto_rico <- as_Spatial(puerto_rico$geometry)
crs(puerto_rico) <- "+proj=longlat +datum=WGS84 +no_defs"

assembly_island_list <- c(list(sao_tome),list(comoros),list(puerto_rico),list(principe))
names(assembly_island_list) <- c("Sao_Tome", "Comoros", "Puerto Rico","Principe")


write_rds(assembly_island_list, file = "Outputs/assembly_islands.rds")
