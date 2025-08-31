\docType{data}
\name{ACT_Covariates}
\alias{ACT_Covariates}
\title{Spatial layers of coarse vegetation types, and distance to a road, trail, or track (metres) for the Australian Capital Territory}
\format{
Two spatial raster layers of a 100 m x 100 m resolution.
\describe{
\item{ForestNGrass}{Coarse vegetation types: forest, grassland and others. This coarse version is derived from raw vegetation types such that forest is a combination of exotic forest, open forest, planation, planting, and woodland; grassland is formed by exotic grassland, grassland, and grassland or secondary woodland. The detail is given in the vignette.}
\item{dist2road_track_trail}{Distance to a road, trail, or track in metres. The distance is calculated using based on the layer of roads, tracks and trails (see Source). }
}
}
\source{
{See inst/extdata/ACT_Covariates.tif in LICENSE.md.}
}
\examples{
ff <- system.file("extdata/ACT_Covariates.tif",package="DSTSOM")
covars <- terra::rast(ff)
}
\description{
A Geotiff (.tif) raster provides two spatial layers of a 100 m x 100 m resolution: coarse vegeation types, and distance to a road, trail or track for the Australian Capital Territory (ACT).  
}
\keyword{datasets}
