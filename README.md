This repository hosts the main code for generating the results for “Land Use and Land Cover Change Significantly Shifted Local Water Availability over past Three Decades”.

There are three folders under this repository:<br>
1.Input datasets in releases.<br>
(1) prep_data.nc for global precipitation derived from MSWEP dataset from 1992 to 2020. <br>
(2) evap_data.nc for global evapotranspiration derived from GLEAM dataset from 1992 to 2020. <br>
(3) LULC_data.nc for land cover types and percentage derived from CCI-LC dataset from 1992 to 2020. <br>
(4) LULCC_class.nc for global major land use and land cover change classes. <br>
(5) WA_LULCC_change.nc for changes in water availability that only caused by LULCC impacts (Data for Fig.4a). <br>
(6) WA_all_change.nc for changes in water availability that caused by total impacts (Data for Fig.4b). <br>
(7) LULCC_weight.nc for the weight of LULCC factor in water availability variations (Data for Fig.4c). <br>
2.LULCC_class.py presents the procedure for the definition of LULCC classes based on CCI satellite dataset. <br>
3.Improve_approximation_method.py presents the space-for-time approximation method for isolating the LULCC impacts on precipitation and evapotranspiration.
