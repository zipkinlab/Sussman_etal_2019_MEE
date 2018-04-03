# A comparative analysis of methods to identify waterbird hotspots (in the Great Lakes)

*Sussman A.L., Gardner B., Adams E.M., Salas L., Kenow K.P., Luukkonen D.R., Monfils M.J., Mueller W.P., Williams K.A., Leduc-Lapierre M., and Zipkin E.F. In review. A comparative analysis of methods to identify species hotspots. Methods in Ecology and Evolution.*

## Abstract
Detection of hotspots is a commonly used method in ecology and conservation to identify areas of high biodiversity or conservation concern. However, delineating and mapping hotspots is subjective and various approaches can lead to different conclusions with regard to the classification of particular areas as hotspots, complicating long-term conservation planning and implementation efforts. We present the results of a comparative analysis of recent approaches for identifying waterbird hotspots. We examined the literature and selected four common measures of identifying persistent areas of high use: two non-parametric, spatial approaches (kernel density estimation and Getis-Ord Gi*) and two parametric, non-spatial approaches (hotspot persistence and hotspots conditional on presence). We applied each of the methods to aerial-survey waterbird count data collected in the Great Lakes from 2012-2014 using a 5 km2 grid. For each approach, we identified areas of high use for seven species/species groups and then compared the results across all methods. Our results indicate that formal hotspot analysis frameworks do not always lead to the same conclusions. The two spatial methods yielded the most similar results across all species analyzed. Yet, we found that the spatial models can differ substantially from the non-spatial models, which were not consistently similar to one another. The hotspot persistence approach differed most significantly from the other methods, but is the only method to explicitly account for temporal variation. We recommend considering the ecological question and scale of any conservation or management activities prior to designing survey methodologies. Deciding the appropriate definition and scale for analysis is also critical for interpretation of hotspot analysis results. Combining methods using an integrative approach (i.e., inclusion of both spatial and non-spatial methods), either within a single analysis or post-hoc, could lead to greater consistency in the identification of waterbird hotspots.

## Data
* alldata_attributed_v04122017.csv - raw aerial visual observations
* data_birds.csv - standardized effort-corrected counts summarized by grid cell (as exported from PostGreSQL)
* data_birds_persistence.csv - standardized effort-corrected counts summarized by day and grid cell (as exported from PostGreSQL)
* data_blocks.csv - unique list of surveyed grid cells (as exported from PostGreSQL)
* ALLSP_poly.zip - zipped shapefile of surveyed grid cells with standardized effort-corrected counts (used for spatial models)

The raw aerial visual observations and transect shapefiles are publicly available through the Midwest Avian Data Center (MWADC), a regional node of the Avian Knowledge Network hosted by Point Blue Conservation Science: https://data.pointblue.org/partners/mwadc/. The observation data were imported into the opensource relational database management system PostgreSQL v9.5.0, with the PostGIS extension v2.2.1, using GDAL ogr2ogr. Shapefiles of survey transects were downloaded from the MWADC and were also imported into PostgreSQL. Data were QA/QC'd and standardized to account for differences in survey methods in PostgreSQL. Using the RPostgreSQL library (Conway et al. 2017), the data were called directly from PostgreSQL within R (RStudio v1.0.136); for archival purposes, the data were exported from PostGreSQL to csv format. Three hotspot models were coded in R, and one model was conducted in ArcGIS v10.3.1.

## Code
* comparative_hotspot_analysis.R - R code to perform three of the four hotspot models (Getis-Ord Gi*, hotspot persistence, and hotspots conditional on presence) used in the comparative analysis; the fourth model (kernel density estimation) was done directly in ArcGIS using the kernel density tool found in the Spatial Analyst toolbox. The models were conducted on the standardized effort-corrected count data. The comparative analysis is also included in this file: a comparison was performed on the results of all four models. Each model's resultant values were first ranked and then a Pearson's correlation was run on the ranked values to determine consistency across the different hotspot models.

## Results (for one example species group: all-species-combined)
**Hotspot models**
* ALLSP_kde.csv - model results for the kernel density estimation approach (as exported from ArcGIS)
* ALLSP_gstat.csv - model results for the Getis-Ord Gi* approach
* ALLSP_persistence.csv - model results for the hotspot persistence approach (for each unique sampling event)
* ALLSP_persistence_totals.csv - hotspot persistence: overall persistence calculated for each block (as exported from PostGreSQL)
* ALLSP_conditional.csv - model results for the hotspots conditional on presence approach

**Comparative analysis**
* ALLSP_vals.csv - model results and grid cells in 75th percentile used in correlations and for mapping
* ALLSP_corr.csv - correlation results

**Species maps**
* pdf file

## Support
This research was supported by the U.S. Fish and Wildlife Service Great Lakes Fish and
Wildlife Restoration Act grants program and was part of a large collaborative effort to monitor and map avian resources
in the Great Lakes. See the Great Lakes Commission project page for more information: https://www.glc.org/work/avian-resources.
