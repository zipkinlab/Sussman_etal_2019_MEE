# [A comparative analysis of common methods to identify waterbird hotspots](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13209)

### Allison L. Sussman, Beth Gardner, Evan M. Adams, Leo Salas, Kevin P. Kenow, David R. Luukkonen, Michael J. Monfils, William P. Mueller, Kate A. Williams, Michele Leduc-Lapierre, and Elise F. Zipkin

### Methods in Ecology and Evolution *In press*

### Code/Data DOI: https://doi.org/10.5061/dryad.rs776p3

### Please contact the first author for questions about the code or data: Allison L. Sussman (sussman7@msu.edu)
__________________________________________________________________________________________________________________________________________
## Abstract
Hotspot analysis is a commonly used method in ecology and conservation to identify areas of high biodiversity or conservation concern. However, delineating and mapping hotspots is subjective and various approaches can lead to different conclusions with regard to the classification of particular areas as hotspots, complicating long-term conservation planning. We present a comparative analysis of recent approaches for identifying waterbird hotspots, with the goal of developing insights about the appropriate use of these methods. We selected four commonly used measures to identify persistent areas of high use: kernel density estimation, Getis-Ord G~i~&ast;, hotspot persistence, and hotspots conditional on presence, which represent the range of quantitative hotspot estimation approaches used in waterbird analyses. We applied each of the methods to aerial survey waterbird count data collected in the Great Lakes from 2012-2014. For each approach, we identified areas of high use for seven species/species groups and then compared the results across all methods and with mean effort-corrected counts. Our results indicate that formal hotspot analysis frameworks do not always lead to the same conclusions. The kernel density and Getis-Ord G~i~&ast; methods yielded the most similar results across all species analyzed and were generally correlated with mean effort-corrected count data. We found that these two models can differ substantially from the hotspot persistence and hotspots conditional on presence estimation approaches, which were not consistently similar to one another. The hotspot persistence approach differed most significantly from the other methods but is the only method to explicitly account for temporal variation. We recommend considering the ecological question and scale of conservation or management activities prior to designing survey methodologies. Deciding the appropriate definition and scale for analysis is critical for interpretation of hotspot analysis results as is inclusion of important covariates. Combining hotspot analysis methods using an integrative approach, either within a single analysis or post-hoc, could lead to greater consistency in the identification of waterbird hotspots.

## Data
* [alldata_attributed_v04122017.csv](https://github.com/zipkinlab/Sussman_etal_MEE/blob/master/alldata_attributed_v04122017.csv) - raw aerial visual observations
* [data_birds.csv](https://github.com/zipkinlab/Sussman_etal_MEE/blob/master/data_birds.csv) - standardized effort-corrected counts summarized by grid cell
* [data_birds_persistence.csv](https://github.com/zipkinlab/Sussman_etal_MEE/blob/master/data_birds_persistence.csv) - standardized effort-corrected counts summarized by day and grid cell
* [data_blocks.csv](https://github.com/zipkinlab/Sussman_etal_MEE/blob/master/data_blocks.csv) - unique list of surveyed grid cells
* [ALLSP_poly.zip](https://github.com/zipkinlab/Sussman_etal_MEE/blob/master/ALLSP_poly.zip) - zipped shapefile of surveyed grid cells with standardized effort-corrected counts (used for spatial models)

The raw aerial visual observations and transect shapefiles are publicly available through the Midwest Avian Data Center (MWADC), a regional node of the Avian Knowledge Network hosted by Point Blue Conservation Science: https://data.pointblue.org/partners/mwadc/. The observation data were imported into the opensource relational database management system PostgreSQL v9.5.0, with the PostGIS extension v2.2.1, using GDAL ogr2ogr. Shapefiles of survey transects were downloaded from the MWADC and were also imported into PostgreSQL. Data were QA/QC'd and standardized to account for differences in survey methods in PostgreSQL. Using the RPostgreSQL library (Conway et al. 2017), the data were called directly from PostgreSQL within R (RStudio v1.0.136); for archival purposes, the data were exported from PostGreSQL to csv format. Three hotspot models were coded in R, and one model was conducted in ArcGIS v10.3.1.

## Code
* [comparative_hotspot_analysis.R](https://github.com/zipkinlab/Sussman_etal_MEE/blob/master/comparative_hotspot_analysis.R) - R code to perform three of the four hotspot models (Getis-Ord G~i~&ast, hotspot persistence, and hotspots conditional on presence) used in the comparative analysis; the fourth model (kernel density estimation) was done directly in ArcGIS using the kernel density tool found in the Spatial Analyst toolbox. The models were conducted on the standardized effort-corrected count data. The comparative analysis (on the results of all four models) is also included in this file. Each model's resultant values were first ranked and then a Pearson's correlation was run on the ranked values to determine consistency across the different hotspot models.

## Results (for one example species group: all-species-combined)
##### Hotspot models
* [ALLSP_kde.csv](https://github.com/zipkinlab/Sussman_etal_MEE/blob/master/ALLSP_kde.csv) - results for the kernel density estimation approach
* [ALLSP_gstat.csv](https://github.com/zipkinlab/Sussman_etal_MEE/blob/master/ALLSP_gstat.csv) - results for the Getis-Ord G~i~&ast approach
* [ALLSP_persistence.csv](https://github.com/zipkinlab/Sussman_etal_MEE/blob/master/ALLSP_persistence.csv) - results for the hotspot persistence approach (for each unique sampling event)
    + [ALLSP_persistence_totals.csv](https://github.com/zipkinlab/Sussman_etal_MEE/blob/master/ALLSP_persistence_totals.csv) - overall persistence calculated for each block
* [ALLSP_conditional.csv](https://github.com/zipkinlab/Sussman_etal_MEE/blob/master/ALLSP_conditional.csv) - results for the hotspots conditional on presence approach

##### Comparative analysis
* [ALLSP_vals.csv](https://github.com/zipkinlab/Sussman_etal_MEE/blob/master/[ALLSP_vals.csv) - model results and grid cells in 75th percentile used in correlations and for mapping
* [ALLSP_corr.csv](https://github.com/zipkinlab/Sussman_etal_MEE/blob/master/ALLSP_corr.csv) - correlation results

##### Species maps
* [ALLSP_maps.pdf](https://github.com/zipkinlab/Sussman_etal_MEE/blob/master/ALLSP_maps.pdf) - example species results maps; includes mean effort-corrected counts, potential hotspots shown on the raw scale for each model, and potential hotspots above the 75th percentile for each model.

## Support
This research was supported by the U.S. Fish and Wildlife Service Great Lakes Fish and
Wildlife Restoration Act grants program and was part of a large collaborative effort to monitor and map avian resources
in the Great Lakes. See the Great Lakes Commission project page for more information: https://www.glc.org/work/avian-resources.
