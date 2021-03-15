<!-- BEGIN COMMENT -->

 [<< Previous Chapter](S5P-WRFChem_UG_ch02_running.md)- [Home](README.md) - [Next Chapter >>](S5P-WRFChem_UG_ch04_analysis.md)

<!-- END COMMENT -->

# 3. Product

## 3.1 Product types

Currently, we have two data products available: one using the WRF-Chem NO2 profile and another using the TM5 profile. Both of them are at the native TROPOMI pixel resolution. Usually, users should care about the first type. If the user want to check some specific improvements, the second one would be useful as the reference.

## 3.2 Version Scheme

The S5P-WRFChem version numbering combines the TROPOMI version number with a type letter. Here's an example:

```
Official product:
	S5P_L2__NO2____20190725T042743_20190725T060912_09219_01_020100_20190911T091424.nc
WRF-Chem product:
	S5P_CHEM_L2__NO2____20190725T042743_20190725T060912_09219_01_020100_20190911T091424.nc
TM5 product:
	S5P_TM5_L2__NO2____20190725T042743_20190725T060912_09219_01_020100_20190911T091424.nc
```

## 3.3 File format

The generated outputs use the same format as the official NetCDF files.

The `/S5P` group contains loaded official variables as DataArrays, while the `/CHEM` or `/TM5` group includes the re-retrieved variables.

*Tools for reading NetCDF file:*

- GUI: [Panoply](https://www.giss.nasa.gov/tools/panoply/)
- Python Script: [xarray](http://xarray.pydata.org/en/stable/) and [netcdf4](https://unidata.github.io/netcdf4-python/)

## 3.4 Variables

The Table below list the attributes of all variables found in the S5P-WRFChem files. Product indicates whether the variable is directly copied from the TROPOMI product (`S5P`), or calculated from S5P-WRFChem (`CHEM/TM5`).

| Varname                                        | Product  | Units                                                        | Description                                                  |
| ---------------------------------------------- | -------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| time                                           | S5P      | days since \<yyyy-mm-dd\>                                    | time using proleptic gregorian calendar                      |
| latitude                                       | S5P      | degrees_north                                                | pixel center latitude                                        |
| longitude                                      | S5P      | degrees_east                                                 | pixel center longitude                                       |
| air_mass_factor_clear                          | S5P      | 1                                                            | Air mass factor for the cloud-free part of the scene         |
| air_mass_factor_cloudy                         | S5P      | 1                                                            | Air mass factor for the cloud-covered part of the scene      |
| air_mass_factor_stratosphere                   | S5P      | 1                                                            | Stratospheric air mass factor                                |
| air_mass_factor_troposphere                    | S5P      | 1                                                            | Tropospheric air mass factor                                 |
| amf_geo                                        | S5P      | 1                                                            | geometric air mass factor<br />amf_geo = 1.0 / cos(solar_zenith_angle) + 1.0 / cos(viewing_zenith_angle) |
| assembled_lat_bounds                           | S5P      | degrees_north                                                | assembled_latitude_bounds calculated by Satpy                |
| assembled_lon_bounds                           | S5P      | degrees_east                                                 | assembled_longitude_bounds calculated by Satpy               |
| cloud_albedo_crb                               | S5P      | 1                                                            | Cloud albedo in the cloud product                            |
| cloud_fraction_crb_nitrogendioxide_window      | S5P      | 1                                                            | Cloud fraction at 440 nm for NO2 retrieval                   |
| cloud_pressure_crb                             | S5P      | Pa                                                           | Cloud optical centroid pressure                              |
| cloud_radiance_fraction_nitrogendioxide_window | S5P      | 1                                                            | Cloud radiance fraction at 440 nm for NO2 retrieval          |
| latitude_bounds                                | S5P      | degrees_north                                                | bounds of latitude                                           |
| longitude_bounds                               | S5P      | degrees_east                                                 | bounds of longitude                                          |
| nitrogendioxide_ghost_column                   | S5P      | mol m-2                                                      | Ghost column NO2: modelled NO2 column below the cloud top    |
| nitrogendioxide_slant_column_density           | S5P      | mol m-2                                                      | Stratospheric vertical column of nitrogen dioxide, derived from the TM5-MP vertical profiles |
| nitrogendioxide_tropospheric_column            | S5P      | mol m-2                                                      | Tropospheric vertical column of nitrogen dioxide             |
| no2_scd_flag                                   | S5P      | 1                                                            | NO2 SCD flag in view of saturation and outliers<br /><br />-1 = no SCD due to saturation limit exceeded, i.e. pqf = 54 ;\n 0 = SCD with Delta <  3.3e-5  &  no error reported ; \n 1 = SCD with Delta <  3.3e-5  &  error reported: pqf=55 (max. outliers exceeded) ;\n 2 = SCD with Delta <  3.3e-5  &  other error reported, e.g. pqf=41 (generic range error) ;\n 3 = SCD with Delta >= 3.3e-5  &  no error reported ;\n 4 = SCD with Delta >= 3.3e-5  &  error reported: pqf=55 ; \n 5 = SCD with Delta >= 3.3e-5  &  other error reported, e.g. pqf=41 ;\nFillValue = no SCD due other error |
| solar_azimuth_angle                            | S5P      | degree<br />clockwise from the North (East = 90, South = 180, West = 270) | Solar azimuth angle at the ground pixel location on the reference ellipsoid. |
| solar_zenith_angle                             | S5P      | degree<br />measured away from the vertical                  | Solar zenith angle at the ground pixel location on the reference ellipsoid. |
| surface_albedo_nitrogendioxide_window          | S5P      | 1                                                            | Surface albedo in the NO2 fit window                         |
| surface_pressure                               | S5P      | Pa                                                           | Surface pressure                                             |
| time_utc                                       | S5P      | 1                                                            | Time of observation as ISO 8601 date-time string             |
| tm5_constant_a                                 | S5P      | Pa                                                           | TM5 hybrid A coefficient at upper and lower interface levels |
| tm5_constant_b                                 | S5P      | Pa                                                           | TM5 hybrid B coefficient at upper and lower interface levels |
| tm5_tropopause_layer_index                     | S5P      | 1                                                            | TM5 layer index of the highest layer in the tropopause       |
| viewing_azimuth_angle                          | S5P      | degree<br />measured clockwise from the North (East = 90, South = 180, West = 270) | Satellite azimuth angle at the ground pixel location on the reference ellipsoid. |
| viewing_zenith_angle                           | S5P      | degree<br />measured away from the vertical                  | Zenith angle of the satellite at the ground pixel location on the reference ellipsoid. |
|                                                |          |                                                              |                                                              |
| amfTrop                                        | CHEM/TM5 | 1                                                            | Tropospheric air mass factor for total tropospheric NO2 calculated by WRF-Chem |
| amfTropVis                                     | CHEM/TM5 | 1                                                            | Tropospheric air mass factor for visible tropospheric NO2 calculated by WRF-Chem |
| swClr                                          | CHEM/TM5 | 1                                                            | Clear-sky NO2 box tropospheric air mass factor based on LUT  |
| swCld                                          | CHEM/TM5 | 1                                                            | Cloudy-sky NO2 box Tropospheric air mass factor based on LUT |
| avgKernel                                      | CHEM/TM5 | 1                                                            | Averaging kernels computed for the weighted average of cloudy and clear conditions |
| no2apriori                                     | CHEM/TM5 | parts-per-part                                               | NO2 a priori profile used for each pixel. Pressure levels given in plevels |
| no2Trop                                        | CHEM/TM5 | mol m-2                                                      | Tropospheric NO2 VCD including estimated ghost column (scdTrop/amfTrop) |
| no2TropVis                                     | CHEM/TM5 | mol m-2                                                      | Tropospheric NO2 VCD with only visible NO2 included (scdTrop/amfTropVis) |
| scdClr                                         | CHEM/TM5 | molec. cm-2                                                  | Clear-sky NO2 slant column density (integration of swClr from the surface to tropopause) |
| scdCld                                         | CHEM/TM5 | molec. cm-2                                                  | Cloudy-sky NO2 slant column density (integration of swCld from the cloud to tropopause) |
| vcdGnd                                         | CHEM/TM5 | molec. cm-2                                                  | NO2 vertical column density (integration of no2apriori from the surface to tropopause) |
| tropopause_pressure                            | CHEM/TM5 | hPa                                                          | TM5 tropopause pressure                                      |
| plevels                                        | CHEM/TM5 | hPa                                                          | TM5 Pressure levels                                          |
| noapriori                                      | CHEM/TM5 | parts-per-part                                               | NO a priori profile used for each pixel. Pressure levels given in plevels |
| vcdGnd_no                                      | CHEM/TM5 | molec. cm-2                                                  | NO vertical column density (integration of noapriori from the surface to tropopause) |

## 3.5 Pixel Filtering

As we don't generate new flags for our products, users should stick with the flags in the `S5P` group.

Two flags are available after version 2: 1) *qa_value* (f<sub>QA</sub>) and 2) *no2_scd_flag*.

Brief explanations are copied here:

#### qa_value (f<sub>QA</sub>)

- 0.75 ≤ f<sub>QA</sub> ≤ 1.00

  The ground pixel is recommended for all applications, including column comparisons, visualisation, trends, monthly/seasonal averages. The data is restricted to cloud-free observations (cloud radiance fraction < 0.5), and snow-ice free observations.

- 0.50 ≤ f<sub>QA</sub> < 0.75

  The ground pixel is recommended for use in data assimilation and comparisons against models or vertical profile observations, given that the averaging kernel is used to specify the sensitivity profile in cloudy situations; this includes good quality retrievals over clouds and snow/ice.

- 0 < f<sub>QA</sub> < 0.50

  f<sub>QA</sub> = 0 The ground pixel is not recommended for use due to serious retrieval issues. A processing error occurred so that the ground pixel cannot be used at all, or the solar zenith angle exceeds the limit set in the data assimilation

#### no2_scd_flag

Only SCD values with flag ‘0’ are in general reliable, though values with flag ‘3’ may still be used with care. But keep here in mind the limit set on the number of allowed outliers for L1B-v1.0 data: we would like to set that at 5 but use for this special processing 15 (first table). Hence it would for these cases be good idea to also take note of the value given in number_of_outliers: if that is above 5, take care, if that is above 10 be even more careful …

SCD values with flags ‘1’ and ‘4’ are highly unreliable as there were too many outliers, so that the DOAS fit is not redone (i.e. the given values are the result of the DOAS fit with spikes in the reflectance!).

SCD values with flag ‘2’ (and ‘5’) might be useful, that is not possible to say beforehand as we do not know where exactly the pqf=41 comes from (this could be related e.g. to reflectance values > 1).

Detailed information is available in the TROPOMI [ATBD document](https://sentinel.esa.int/web/sentinel/technical-guides/sentinel-5p/products-algorithms).

## 3.6 Column Densities

Like [BEHR](https://github.com/CohenBerkeleyLab/BEHR-core), the field *no2Trop* contains the to-ground column, which includes the estimated ghost column below the cloudy part of the pixel. The field *no2TropVis* does not include the ghost column; only the above-cloud NO<sub>2</sub> is included for the cloudy component of the pixel.

- Most users should use the *no2Trop* field.
- Users interested in cloud slicing approaches should use the *no2TropVis* field.

Details are shown in the [BEHR User's Guide](https://github.com/CohenBerkeleyLab/BEHR-core/tree/develop/Documentation/BEHR%20User%20Guide).

<!-- BEGIN COMMENT -->

 [<< Previous Chapter](S5P-WRFChem_UG_ch02_running.md)- [Home](README.md) - [Next Chapter >>](S5P-WRFChem_UG_ch04_analysis.md)<br>
S5P-WRFChem User's Guide (c) 2021<br>

<!-- END COMMENT -->