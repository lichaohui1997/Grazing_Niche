# Grazing_Niche


Read me for the code:

Code to reproduce the analysis and figures for the article: CH. Li, M. Kotz, P. Pradhan, XD. Wu, YC. Hu, Z. Li, GQ. Chen. Global Climate Niche for Grazing Implies Shifting Suitability    

Data Primary data are not provided in the repository due to storage capacity. Full download of these data requires storage capacity of 246GB. All primary data used in this study are openly accessible at the following sites:  

CMIP (Coupled Model Intercomparison Project) can be accessed at https://pcmdi.llnl.gov/CMIP5/.  

ESFG data (Earth System Grid Federation) can be accessed at https://esgf.llnl.gov CHELSA V2 (Climatologies at high resolution for the earth’s land surface areas) can be accessed at https://chelsa-climate.org.  

Hyde (History database of the Global Environment) can be accessed at https://www.pbl.nl/en/hyde-history-database-of-the-global-environment.  

Gridded Livestock of the World (GLW 3) database can be accessed at https://livestockdata.org/contributors/food-and-agriculture-organization-united-nations/gridded-livestock-world-glw3.  

Gridded Population of the World Version 4 (GPWv4) from SEDAC (Socioeconomic Data and Applications Center) can be accessed at https://sedac.ciesin.columbia.edu/data/collection/gpw-v4.  

The MOD17A3HGF Version 6 dataset can be accessed at https://lpdaac.usgs.gov/products/mod17a3hgfv006/.  

Percent Tree Coverage (PTC) Global version can be accessed at https://globalmaps.github.io/ptc.html.  

EARTH Env data can be accessed at https://www.earthenv.org. Land-Use Harmonization (LUH2) can be accessed at https://luh.umd.edu.  

Global Aboveground and Belowground Biomass Carbon Density Maps can be accessed at https://daac.ornl.gov/VEGETATION/guides/Global_Maps_C_Density_2010.html.  

SOC data is from literature (https://pnas.org/doi/full/10.1073/pnas.1706103114.  


To replicate the results, please run the 48 code scripts in subsequent orders from run1 series to run24 series. 

Script run1_historical_simulation.m downloads and process the historical climate data for the variable 'pr', 'tas', 'hurs', 'sfcWind' from CMIP, and land use data from HYDE database. 

Script run1_historical_simulation_regional.m conducts region-specific reanalysis using historical data, run1_sensitivity_timeframe.m tests the sensitivity of the historical timeframe chosen.

Script run2_2newclimate.m and run3_5historical_climate.m add more climate data to the ensemble and conduct niche analysis using historical climate and grazing data.

Script run4_modern_climateniche.m downloads and processes modern datasets, including downscaled CHELSA datasets

Script run5_livestock.m process data on modern livestock distribution data 

Script run5_montecarlo.m conducts monte-carlo analysis

Script run6_future.m downloads and processes future climate datasets and land use datasets to estimate future niches.
run6_future_regional.m conducts continental re-analysis with a similar procedure to run6_future.m.

Script run7_analysis.m performs analysis on above-ground biomass, below-ground biomass, and overgrazing analysis

Script run8_analysis_future.m performs further niche analysis in future climate scenarios

Script run8_analysis_future_regional.m performs regional niche analysis for future climate scenarios

Script run9_figures_all.m creates figures 

Script run10_modernniche.m performs analysis using modern livestock distribution data to determine the niche

Script run11_latlon_future.m and run12_latlon_modern.m performs longitudinal and latitudinal analysis 

Script run13_figures_allnew.m creates more figures

Script run14_presentfuture_climate_compare.m compares the differences between the temperature between 2100 and present

Script run15_latitude.m calculates the latitude dispersion of grassland for present and future

Script run16_impacted_poplivestock_2100.m calculates the turnover impact on population and livestock figure

Script run16_impacted_poplivestock_2100_sensitivity.m calculates the uncertainty for the impacted population using the widest and thinnest thresholds

Script run17_sensitiviti_all_cattle.m, run17_sensitiviti_all_sheep.m, and run17_sensitivity_all_goats.m calculates the sensitivity analysis for different livestock species

Script run17_sensitivity_all_turnover.m calculates the turnover rate of future and present niche

Script run18_continental_presentniche.m conducts continental re-analysis, as for all run18 series scripts.

Script run19_continental_futureniche.m performs reanalysis for continental turnover.

Script run20 calculates the reanalysis for using alternate CMIP dataset to replace CHELSA dataset

Script run21_newniche_MRIdata.m conducts sensitivity analysis for using alternative climate datasets to calculate niche

Script run23_presentfuture_climate_compare_MRI.m compares the differences between the temperature between 2100 and present using the same set of MRI data

Script run22_turnover_MRIdata.m uses alternative climate dataset to calculate the turnover

Script run24_species_specific_heatmap.m produces species-specific heatmaps

Script run25_treemap_2100.R creates heatmaps in Fig. 3


Key scripts:
run1-run3 gathers and processes data of historical land use and historical climate.
run13_figures_allnew has the scatter plots that are in Figure 1.
run10 calculates the livestock niche, the cattle niche, sheep niche, and goats niche
run17_sensitiviti_all_turnover calculates the turnover of average, wide, and thinnest thresholds. This results in figure2
run16_impacted_poplivestock_2100 and run16_impacted_poplivestock_2100_uncertainty calculate the results for figure3, which results in the impacted population with uncertainty ranges

Additional key sensitivity analysis scripts:
run18 calculatest the continental niche
run20 calculates the reanalysis for using alternate CMIP dataset to replace CHELSA dataset
run24 produces species-specific heatmaps
run25 creats treemap in Figure 3

Dependencies required   
For Matlab, our scripts require the dependencies require the climate toolbox (https://zenodo.org/badge/latestdoi/171331090), land mask, and addcolorplus package. All scripts were run in version Matlab_R2021b environment. For any further questions please feel free to contact: Li.Chaohui@pik-potsdam.de  
