'''
Regrid WRF-Chem IRR outputs to S5P pixels

UPDATE:
    Xin Zhang:
       01/08/2022: Basic
'''

import os
import logging
import warnings
import pandas as pd
import xarray as xr
from s5p_functions import Config
from distributed import Client, LocalCluster
from s5p_utils import validation, load_s5p, load_wrf, load_irr, regrid_irr, save_irr

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
            s5p, _, _, mean_t = load_s5p(curr_date, cfg['s5p_nc_dir'])
            output_file = os.path.join(cfg['output_data_dir'], 'S5P_IRR'+s5p.attrs['s5p_filename'][8:])

            # skip if the file exists
            if os.path.isfile(output_file) and not overwrite:
                logging.info('File already exists, skipping')
                continue

            # load wrf data
            wrf_file, wrf = load_wrf(mean_t, cfg['wrf_nc_dir'])
            # load irr data
            irr = load_irr(mean_t, cfg['wrf_nc_dir'])

            # regrid IRR data to pixels
            interp_da_list = []

            for index,t in enumerate(irr.coords['Time']):
                logging.info(f"Processing IRR data at {t.dt.strftime('%Y-%m-%d %H:%M').values}")
                interp_da = regrid_irr(irr.sel(Time=t), wrf, s5p)
                interp_da_list.append(interp_da.load())
                del interp_da

            # concat data into one DataArray by "Time" dimension
            ds_irr = xr.concat(interp_da_list, dim='Time')

            # save data
            save_irr(s5p, ds_irr, output_file)


if __name__ == '__main__':
    main()
