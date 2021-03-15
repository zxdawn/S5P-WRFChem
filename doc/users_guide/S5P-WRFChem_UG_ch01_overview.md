<!-- BEGIN COMMENT -->

[Home](README.md) - [Next Chapter >>](S5P-WRFChem_UG_ch02_running.md)

<!-- END COMMENT -->

# 1. Overview

## 1.1 Introduction

S5P-WRFChem was initiated by **Xin Zhang**, a Ph. D. candidate in NUIST. He showed that using high-resolution WRF-Chem simulations could improve the NO<sub>2</sub> retrieval (especially lightning NO<sub>x</sub>) based on the TROPOMI NO<sub>2</sub> Standard Product. You can check his papers in the Literatures section below or [his blog](https://dreambooker.site/) for more information.
If you want to contribute to S5P-WRFChem, add your name to the [AUTHOR list](https://github.com/zxdawn/S5P-WRFChem/blob/master/AUTHORS.md). The bold name is the most recent contact if you need help.

## 1.2 Literatures

The following are papers using S5P-WRFChem algorithm:

- Coming soon ...

## 1.3 Features

- Re-calculation of tropospheric NO<sub>2</sub> Air Mass Factors (AMFs) using WRF-Chem simulations
- Re-retrieval of tropospheric NO<sub>2</sub> column density
- Retrieval of lightning NO<sub>x</sub>

## 1.4 Improvements (to-do)

- Manual pressure levels

  This could utilize the high-resolution WRF-Chem outputs better, especially for convections.

- Fixed lon/lat grids

  The fixed lon/lat grid can be useful for long term reanalysis. This can be achieved by [xESMF](https://xesmf.readthedocs.io/en/latest/) or [Pyresample](https://pyresample.readthedocs.io/en/latest/) packages.

- New datasets for snow/ice flag

  TROPOMI version <1.4 used the NICE database and switchs to the ECMWF snow/ice data. If there're more accurate data, it's better to add the support.

<!-- BEGIN COMMENT -->

[Home](README.md) - [Next Chapter >>](S5P-WRFChem_UG_ch02_running.md)<br>
S5P-WRFChem User's Guide (c) 2021<br>

<!-- END COMMENT -->