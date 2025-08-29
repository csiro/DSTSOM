\name{ACT_FloristicSurveys}
\alias{ACT_FloristicSurveys}
\title{A dataset of floristics data, weed control, and covariates for the Australian Capital Territory (ACT)}
\format{
A data frame of 1299 observations containing information on
  the following 14 variables:
\describe{
  \item{year}{Monitoring year}
  \item{easting}{The easting (in metres) of a subplot}
  \item{northing}{The northing (in metres) of a subplot}
  \item{eragrostis_curvula_0}{The modified Braun-Blanquet score for African
    lovegrass from \code{0} to \code{7}: Score \code{0} represents 0\%
    cover; \code{1} is \eqn{<}{<} 5\% cover and solitary (\eqn{<}{<} 4 individuals); \code{2}
    is \eqn{<}{<} 5\% cover and few (4 to 15 individuals); \code{3} is \eqn{<}{<} 5\% cover
    and few (\eqn{>}{>} 15 individuals); \code{4} is \eqn{>=}{>=} 5 \% cover  \eqn{<}{<} 25 \%;
    \code{5} is \eqn{>=}{>=} 25 \% cover  \eqn{<}{<} 50 \%; \code{6} is \eqn{>=}{>=} 25 \% cover  \eqn{<}{<} 50 \%;
    \code{6} is \eqn{>=}{>=} 50 \% cover  \eqn{<}{<} 75 \%; \code{6} is 75 \% cover
    and greater.}
  \item{nassella_trichotoma_0}{The modified Braun-Blanquet score for serrated tussock.}
  \item{nassella_neesiana_0}{The modified Braun-Blanquet score for Chilean needle grass.}
  \item{hypericum_perforatum_0}{The modified Braun-Blanquet score for  St. John's wort.}
  \item{ALG_daySinceControl}{The number of days since the last control of African lovegrass. This is derived based on the monitoring time information in the floristic data collected by the ACT monitoring team and the control time information in the control data collected by the ACT control team.}
  \item{ST_daySinceControl}{The number of days since the last control of serrated tussock.}
  \item{CNG_daySinceControl}{The number of days since the last control of Chilean needle grass.}
  \item{StJW_daySinceControl}{The number of days since the last control of St. John's wort.}
  \item{veg_structure}{Vegetation types}
  \item{dist2road_track_trail}{Distance to a road, trail, or track in metres. The distance is calculated based on the layer of roads, tracks and trails.}
  \item{ForestNGrass}{Coarse vegetation types: forest, grassland and others. This coarse version is derived from veg_structure such that forest is a combination of exotic forest, open forest, planation, planting, and woodland; grassland is formed by exotic grassland, grassland, and grassland or secondary woodland. The detail is given in the vignette.}
}
}
\source{
{See Data/ACT_FloristicSurveys.RData in LICENSE.md.}
}
\usage{
data(ACT_FloristicSurveys, package = "DSTSOM")
}
\description{
The data table consists of floristics data from the ACT
  monitoring team, weed control metadata from the ACT weed control
  team, and relevant covariates.
}
\keyword{datasets}

