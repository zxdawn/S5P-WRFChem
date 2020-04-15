'''
INPUT:
    WRF-Chem files
    S5p (TROPOPMI L2) files
    Look-up table (LUT)

OUTPUT:
    S5p data retrieved by wrf-chem no2 profiles

STEPS:
    1. Iterate through dates
    2. Read s5p (TROPOMI L2) data by "Satpy"
    3. Read wrf-chem files nearest to s5p swath
    4. Calculate the AMFs
        4.1 Interpolate wrf profiles to pressure (set in settings.txt)
        4.2 Average (area_weighted) profiles to s5p pixels by "xESMF"
        4.3 Use LUT to calculate Box-AMFs:
                bamf_clr and bamf_cld
        4.4 Use Box-AMFs to calculate AMFs:
                amf, amfVis, amfLNO2 and amfLNOx
    5. Calculate no2Trop, no2TropVis, LNO2 and LNOx
    6. Save to NetCDF files

UPDATE:
    Xin Zhang:
       04/12/2020: Basic

'''

import os
import logging
import pandas as pd
from s5p_utils import *
from s5p_plots import *
from distributed import Client, LocalCluster

logging.getLogger('satpy').setLevel(logging.ERROR)
# Choose the following line for info or debugging:
logging.basicConfig(level=logging.INFO)
# logging.basicConfig(level=logging.DEBUG)


def main():
    # read config file and validate it
    cfg = Config('settings.txt')
    logging.info(cfg)
    overwrite = cfg.get('overwrite', 'True') == 'True'
    plot_weights = cfg.get('plot_weights', 'True') == 'True'
    plot_regrid = cfg.get('plot_regrid', 'True') == 'True'
    plot_bamf = cfg.get('plot_bamf', 'True') == 'True'
    plot_interp = cfg.get('plot_interp', 'True') == 'True'
    # pres_levels = np.array([int(x) for x in cfg['pres_levels'].split(',')])
    validation(cfg)

    # set cluster
    cluster = LocalCluster(n_workers=int(cfg['n_workers']),
                           threads_per_worker=int(cfg['threads_per_worker']),
                           processes=True,)
    client = Client(cluster)
    with Client(cluster) as client:
        # generate data range
        req_dates = pd.date_range(start=cfg['start_date'],
                                  end=cfg['end_date'],
                                  freq='D')

        # iterate through dates
        for curr_date in req_dates:
            logging.info(f'Now processing {curr_date.date()} data ...')

            # load s5p data
            s5p, vnames, lut, mean_t = load_s5p(curr_date, cfg['s5p_nc_dir'])

            # create output name based on date
            output_file = os.path.join(cfg['output_data_dir'],
                                       'S5P_CHEM' + s5p.attrs['s5p_filename'][8:])
            if os.path.isfile(output_file) and not overwrite:
                logging.info('File already exists, skipping')
                continue

            # load wrf data
            wrf_file, wrf = load_wrf(mean_t, cfg['wrf_nc_dir'])

            # regrid wrf data to pixels
            regridder, interp_ds = regrid(wrf, s5p)

            # calculate box amf
            bAmfClr, bAmfCld = cal_bamf(s5p, lut, interp_ds)

            # calculate amf, scattering weights and averaging kernel
            ds = cal_amf(s5p, interp_ds, bAmfClr, bAmfCld)

            # save to nc file
            save(s5p, wrf, ds, vnames, output_file)

    # Check these paras at the end of date
    # Plot the weights once
    if plot_weights:
        logging.info('Plot weights_sparse and weights_heatmap')
        plot_weights_sparse(regridder, cfg['output_fig_dir'])
        plot_weights_heatmap(wrf, s5p, regridder, wrf_file, cfg['output_fig_dir'])

    if plot_regrid:
        logging.info('Plot regridded pressure')
        plot_regrid_p(wrf, interp_ds, cfg['output_fig_dir'])

    if plot_interp:
        logging.info('Plot interpolated profile')
        plot_interp_p(interp_ds, cfg['output_fig_dir'])

    # Plot the bamf once
    if plot_bamf:
        logging.info('Plot Box-AMFs')
        plot_bamf_p(bAmfClr, bAmfCld, cfg['output_fig_dir'])


if __name__ == '__main__':
    main()
