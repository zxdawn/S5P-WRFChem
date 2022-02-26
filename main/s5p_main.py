'''
INPUT:
    - WRF-Chem output files
    - S5P (TROPOPMI) L2 product files
    - S5P Look-up table (LUT)

OUTPUT:
    S5P data retrieved by wrf-chem NO2 profiles

STEPS:
    1. Iterate through input dates
    2. Read S5P data by "Satpy" package
    3. Read wrf-chem files nearest to s5p swath
    4. Calculate the AMFs
        4.1 Interpolate wrf profiles to TM5 pressure levels
        4.2 Average (area_weighted) profiles to S5P pixels by "xESMF"
        4.3 Use LUT to calculate Box-AMFs
        4.4 Use Box-AMFs to calculate AMFs
    5. Calculate avgKernel, no2Trop, no2TropVis, and vcdGnd
    6. Save all variables (set in attrs_table.csv) as NetCDF files

UPDATE:
    Xin Zhang:
       04/12/2020: Basic
       03/09/2021: Fix AMF bugs
'''

import os
import logging
import warnings
import pandas as pd
from s5p_utils import *
from s5p_plots import *
from distributed import Client, LocalCluster

logging.getLogger('satpy').setLevel(logging.ERROR)
# Choose the following line for info or debugging:
logging.basicConfig(level=logging.INFO)
# logging.basicConfig(level=logging.DEBUG)

# Disable a few warnings:
warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)
warnings.filterwarnings('ignore', category=UserWarning, append=True)
warnings.filterwarnings('ignore', category=FutureWarning, append=True)


def main():
    # read config file and validate it
    cfg = Config('settings.txt')
    logging.info(cfg)

    overwrite = cfg.get('overwrite', 'True') == 'True'
    tm5_prof = cfg.get('tm5_prof', 'True') == 'True'
    plot_weights = cfg.get('plot_weights', 'True') == 'True'
    plot_comp = cfg.get('plot_comp', 'True') == 'True'
    plot_bamf = cfg.get('plot_bamf', 'True') == 'True'

    validation(cfg)

    # set cluster
    cluster = LocalCluster(n_workers=int(cfg['n_workers']),
                           threads_per_worker=int(cfg['threads_per_worker']),
                           processes=True,)

    with Client(cluster) as client:
        # generate data range
        req_dates = pd.date_range(start=cfg['start_date'],
                                  end=cfg['end_date'],
                                  freq='D')

        # iterate through dates
        for index, curr_date in enumerate(req_dates):
            logging.info(f'Now processing {curr_date.date()} data ...')

            # load s5p data
            s5p, vnames, lut, mean_t = load_s5p(curr_date, cfg['s5p_nc_dir'], tm5_prof)

            # create output name based on date
            if tm5_prof:
                prefix = 'S5P_TM5'
            else:
                prefix = 'S5P_CHEM'
            output_file = os.path.join(cfg['output_data_dir'],
                                       prefix + s5p.attrs['s5p_filename'][8:])

            # skip if the file exists
            if os.path.isfile(output_file) and not overwrite:
                logging.info('File already exists, skipping')
                continue

            # load wrf data
            wrf_file, wrf = load_wrf(mean_t, cfg['wrf_nc_dir'])

            # regrid wrf data to pixels
            regridder, interp_ds = regrid(wrf, s5p, tm5_prof)

            # calculate box amf
            s5p_origin, bAmfClr, bAmfCld, del_lut = cal_bamf(s5p, lut)

            # calculate amf, scattering weights and averaging kernel
            ds = cal_amf(s5p, s5p_origin, interp_ds, bAmfClr, bAmfCld, del_lut)

            # save to nc file
            save(s5p, wrf, ds, vnames, output_file, tm5_prof)

            if index == len(req_dates) -1:
                '''Check these paras at the end of date'''
                # Plot the resample weights
                if plot_weights and regridder:
                    logging.info('Plot weights_sparse and weights_heatmap')
                    plot_weights_sparse(regridder, cfg['output_fig_dir'])
                    plot_weights_heatmap(wrf, s5p, regridder, wrf_file, cfg['output_fig_dir'])

                # compare new retrieval with the official product 
                if plot_comp:
                    logging.info('Plot comparison between S5P official product and CHEM')
                    plot_comp_s5p_chem(s5p, ds, cfg['output_fig_dir'])
        
                # Plot the bamf once
                if plot_bamf:
                    logging.info('Plot Box-AMFs')
                    plot_bamf_p(bAmfClr, bAmfCld, cfg['output_fig_dir'])


if __name__ == '__main__':
    main()

