\docType{data}
\name{CNG_DaysSinceControl}
\alias{CNG_DaysSinceControl}
\title{Spatial layers of days since control for Chilean needle grass in the Australian Capital Territory (ACT)}
\format{
Five spatial layers represent days since the last Chilean needle grass control, referencing July 1st of a given year at a resolution of 100 m x 100 m, from 2019 to 2023.
\describe{
\item{DaySinceCon_Ref20190701}{days since last control including on 1 July 2019}
\item{DaySinceCon_Ref20200701}{days since last control including on 1 July 2020}
\item{DaySinceCon_Ref20210701}{days since last control including on 1 July 2021}
\item{DaySinceCon_Ref20220701}{days since last control including on 1 July 2022}
\item{DaySinceCon_Ref20230701}{days since last control including on 1 July 2023}
}
}
\source{
{See inst/extdata/CNG_DaysSinceControl.tif in LICENSE.md.}
}
\usage{
ff <- system.file("extdata/CNG_DaysSinceControl.tif",package="DSTSOM")
ALGcontrol <- terra::rast(ff)
}
\description{
A Geotiff (.tif) raster provides yearly spatial layers of days since control for Chilean needle grass in a 100 m x 100 m resolution in ACT. The spatial layers are derived based on the monitoring time information in the floristic data collected by the ACT monitoring team and the control time information in the control data collected by the ACT control team.  
}
\keyword{datasets}
