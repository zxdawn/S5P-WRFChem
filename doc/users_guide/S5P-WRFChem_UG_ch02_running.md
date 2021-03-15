<!-- BEGIN COMMENT -->

[<< Previous Chapter](S5P-WRFChem_UG_ch01_overview.md) - [Home](README.md) - [Next Chapter >>](S5P-WRFChem_UG_ch03_product.md)

<!-- END COMMENT -->

# 2. Running S5P-WRFChem

## 2.1 System Recommendation

S5P-WRFChem is a simple retrieval tool written in Python. S5P-WRFChem execution is typically performed on Linux based systems. The hardware configuration of such a system depends on the simulation duration. You can run the retrieval of one day and check whether the disk space is enough.

## 2.2 Python Environment

1. [Clone](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository) or download this repository.

2. Download the latest anaconda package from the [official website](https://www.anaconda.com/) or [TUNA (quicker for China region)](https://mirrors.tuna.tsinghua.edu.cn/anaconda/archive/?C=M&O=A).

3. Install anaconda by running the `.sh` file.

4. Add Anaconda to the system environment following this [tutorial](https://www.pythonlikeyoumeanit.com/Module1_GettingStartedWithPython/Installing_Python.html) or any other online tutorial

5. Check whether the Anaconda environment is set successfully

   ```
   which conda
   ```

6. Changing the path to the root directory `S5P-WRFChem` and create the new environment in the terminal

   ```
   conda env create -f environment.yml
   ```

7. Activate the environment

  ```
  conda activate s5p-wrfchem
  ```

## 2.3 Data Preparation

### TROPOMI

- [SCIHUB](https://scihub.copernicus.eu/)

  This is the link on the official [TROPOMI website](http://www.tropomi.eu/data-products/nitrogen-dioxide). It's convenient to download selected data using **browser** and DownThemAll plugin.

- [GES DISC](https://daac.gsfc.nasa.gov/?)

  The NASA website offers the same version TROPOMI products, but it provides the instruction of downloading data using script.

  This is useful for **long-term analysis and Linux system users**.

- [TEMIS](https://www.temis.nl/airpollution/no2.php)

  The KNMI website is simple to quick view and download the **daily and monthly** TROPOMI product. However, it doesn't support define your own area to limit the file counts.

### WRF-Chem

The [WRF-Chem website](https://www2.acom.ucar.edu/wrf-chem) has provided the detailed User guide and the retrieval algorithm supports wrf output files of all versions.

Note that you need to modify `frames_per_outfile` to `1` in `namelist.input` for the domain which is used for the retrieval. This will generate one file every `history_interval`.

**Put the TROPOMI and WRF-Chem data in the directories set by `setting.txt` in the `main` directory.**

## 2.4 Retrieval

1. Change the directory to the `main` directory and modify the `setting.txt` file. You can just edit the `start_date` and `end_date` at this stage. Feel free to modify other parameters following the instructions.
2. Run `s5p_main.py` in the `main` directory with the new Python Environment: `python s5p_main.py`

<!-- BEGIN COMMENT -->

[<< Previous Chapter](S5P-WRFChem_UG_ch01_overview.md) - [Home](README.md) - [Next Chapter >>](S5P-WRFChem_UG_ch03_product.md)

S5P-WRFChem User's Guide (c) 2021<br>

<!-- END COMMENT -->