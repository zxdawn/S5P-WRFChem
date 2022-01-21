'''
Some functions for s5p_utils.py

UPDATE:
    Xin Zhang:
       04/12/2020: Basic

'''

import os
import scipy
import logging
import numpy as np
import pandas as pd
import xarray as xr
import scipy.interpolate
from datetime import datetime, date
from configparser import SafeConfigParser


class Config(dict):
    def __init__(self, filename):
        ''' Reads in the settings.txt file

        Example 'settings.txt':
        [ALL]
        start_date = 2019-07-25     ; The beginning of s5p data (yyyy-mm-dd)
        end_date = 2019-07-25       ; The end of s5p data (yyyy-mm-dd)
        s5p_nc_dir = s5p_data       ; where s5p data are saved
        wrf_nc_dir = wrf_data       ; where wrfchem data are saved
        output_data_dir = output    ; where generated files are saved
        output_fig_dir = figures    ; where figures are saved
        overwrite = False           ; whether overwrite the same generated files (True/False)
        tm5_prof = False            ; whether read the TM5-MP a-priori profile (True/False)
        n_workers = 4               ; Dask: number of workers to start
        threads_per_worker = 1      ; Dask: number of threads per each worker
        plot_weights = False        ; whether plot the weights
        plot_comp = False           ; whether plot the comparison between S5P and CHEM product
        plot_regrid = False         ; whether plot the regridded surface pressure
        plot_bamf = False           ; whether plot the bamfs
        '''

        curr_dir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))
        config = SafeConfigParser(inline_comment_prefixes=';')
        config.read(filename)
        for key, value in config.items('ALL'):
            if 'dir' in key:
                self[key] = os.path.join(curr_dir, value)
            else:
                self[key] = value


def validate_path(path_in, var_name):
    '''Validate input path'''
    if not isinstance(path_in, str):
        raise ValueError('{} must be a string'.format(var_name))
    elif not os.path.isdir(path_in):
        print(path_in)
        raise ValueError('{} ({}) does not exist'.format(var_name, path_in))


def validate_date(date_in):
    '''Validate input date'''
    try:
        datetime.strptime(date_in, '%Y-%m-%d')
    except ValueError:
        raise ValueError("Incorrect data format, should be YYYY-MM-DD")


def construct_filename(file_date):
    '''Construct the name for output files'''
    if not isinstance(file_date, (datetime, date)):
        raise InputError('file_date must be an instance of '
                         'datetime.date or datetime.datetime')

    return 'S5P_CHEM_{file_date.strftime("%Y%m%d")}.nc'


def interp1d_sp(data, x, xi):
    '''Linear interpolate function'''
    f = scipy.interpolate.interp1d(x, data, fill_value='extrapolate')

    return f(xi)


def xr_interp(input_array,
              input_p, input_p_dimname,
              interp_p, interp_p_dimname):
    '''Interpolate 3D array by another 3D array

    Args:
        input_array:
                the original array
                - dims: input_p_dimname, y, x
        input_p:
                pressure of the original array
                - dims: input_p_dimname, y, x
        interp_p:
                the pressure levels which input_array is interpolated to
                - dims: interp_p_dimname, y, x
        input_p_dimname:
                the name of the vertical dim for input_p
        interp_p_dimname:
                the name of the vertical dim for interp_p

    '''

    logging.debug(' '*8 + f'Interpolating from {input_p_dimname} to {interp_p_dimname}')

    return xr.apply_ufunc(
        interp1d_sp,
        input_array.chunk({
                          input_p_dimname: input_array.sizes[input_p_dimname],
                          'y': input_array.shape[1],
                          'x': input_array.shape[2],
                          }),
        input_p,
        interp_p,
        input_core_dims=[[input_p_dimname], [input_p_dimname], [interp_p_dimname]],
        output_core_dims=[[interp_p_dimname]],
        exclude_dims=set((input_p_dimname,)),
        vectorize=True,
        dask="parallelized",
        output_dtypes=[input_array.dtype],
    )


def interp_to_tm5(regrid_vars, s5p):
    '''Interpolate regridded data to TM5 pressure levels

    References:

        http://xarray.pydata.org/en/stable/examples/apply_ufunc_vectorize_1d.html
        https://github.com/pydata/xarray/issues/3931

    '''

    logging.info(' '*4 + 'Interpolating gridded WRF-Chem data to TM5 pressure levels ...')

    # remove "p" as we don't need to interpolate the pressure data
    varnames = list(regrid_vars.keys())
    varnames.remove('p')

    # combine DataArrays into one DataArray with expanded dimension
    combined = xr.concat(regrid_vars.values(),dim=list(regrid_vars.keys())).rename({'concat_dim':'varname'}).rename('profiles')

    # interpolate data to s5p pressure levels
    interp_da = np.exp(xr_interp(np.log(combined.drop_sel(varname='p')),
                                          np.log(combined.sel(varname='p')), combined.sel(varname='p').dims[0],
                                          np.log(s5p['p'].rolling({'layer': 2}).mean()[1:, ...].load()), 'layer')
                       .transpose('layer',
                                  ...,
                                  transpose_coords=False)
                       )

    # convert to ds for simplier index later
    interp_ds = interp_da.to_dataset(dim='varname')

    logging.info(' '*8 + 'Finish interpolation')

    return interp_ds


def interp_to_layer(profile, pclr, pcld):
    '''Interpolate data to pressure levels including additional layer (surface/cloud)'''

    # calculate the full level pressure which is used by no2
    pclr = pclr.rolling({pclr.dims[0]: 2}).mean()[1:, ...]
    pcld = pcld.rolling({pcld.dims[0]: 2}).mean()[1:, ...]

    return np.exp(xr_interp(np.log(profile),
                            np.log(pclr), pclr.dims[0],
                            np.log(pcld), pcld.dims[0]).transpose(
                            pcld.dims[0],
                            ...,
                            transpose_coords=False))


def concat_p(layers, ptm5):
    '''Concatenate surface pressures, cloud pressures and tm5 pressures'''
    s5p_pcld = xr.concat([ptm5, xr.concat(layers, 'layer')], dim='layer')

    # -- bak ---
    # # https://github.com/pydata/xarray/issues/3954
    # pcld = pcld.expand_dims(layer=[ptm5.sizes['layer']]).drop(['layer'])
    # s5p_pcld = xr.concat([ptm5, pcld], dim='layer')
    # -- bak ---

    # --- pure xarray method ---
    logging.info(' '*6 + 'Sorting pressures ...')
    # get the sorting index
    sort_index = (-1*s5p_pcld.load()).argsort(axis=0)  # (layer, y, x)
    # sort pressures
    s5p_pcld = xr.DataArray(np.take_along_axis(s5p_pcld.values, sort_index, axis=0),
                            dims=['plevel', 'y', 'x'])
    # assign plevel coordinates
    s5p_pcld = s5p_pcld.assign_coords(plevel=range(s5p_pcld.sizes['plevel']))
    # --- pure xarray method ---

    # # --- pandas method (bak) ---
    # # sort pressures in DataFrame
    # # sort_values isn't available for dask DataFrame yet
    # #   https://github.com/dask/dask/issues/958
    # logging.info(' '*6 + 'Sorting pressures slowly...')
    # df = s5p_pcld.to_series()
    # df = df.reset_index() \
    #        .sort_values(['y', 'x', 'tm5_pressure'], ascending=False) \
    #        .set_index(['layer', 'y', 'x'])

    # # Unfortunately, `set_levels` has bug ...
    # # Check issues submitted by Xin:
    # #   https://github.com/pandas-dev/pandas/issues/33420
    # #   https://github.com/pydata/xarray/issues/3957
    # # df.index.set_levels(np.arange(s5p_pcld.sizes['layer']),
    # #                     level='layer',
    # #                     inplace=True)

    # # set new indexes for 'layer'
    # df.index = pd.MultiIndex.from_arrays(
    #                 [np.tile(np.arange(s5p_pcld.sizes['layer']),
    #                          s5p_pcld.sizes['y']*s5p_pcld.sizes['x']),
    #                  df.index.get_level_values(1),
    #                  df.index.get_level_values(2)]) \
    #              .set_names('plevel', level=0)

    # # convert DataFrame to dask DataArray
    # s5p_pcld = df.to_xarray()['tm5_pressure'] \
    #               .chunk({'plevel': s5p_pcld.shape[0],
    #                       'y': s5p_pcld.shape[1],
    #                       'x': s5p_pcld.shape[2]})

    # # we need to drop coords for y and x first
    # #   https://github.com/pydata/xarray/issues/2283
    # s5p_pcld = s5p_pcld.drop(['y', 'x'])#.rename({'plevel': 'layer'})
    # # --- pandas method ---

    return s5p_pcld


def merge_da(profile, p, psfc, ptropo):
    '''Merge DataArrays to Dataset'''
    profile = profile.to_dataset(name='no2')

    return xr.merge([profile, p, psfc.rename('surface_pressure'), ptropo])


def cal_tropo(pclr, itropo):
    '''Calculate tropopause pressure based on index'''

    # get fill_values and overwrite with 0, because index array should be int.
    # Check submitted issue by Xin:
    #   https://github.com/pydata/xarray/issues/3955
    tropo_bool = itropo == itropo._FillValue
    itropo = xr.where(tropo_bool, 0, itropo)

    # isel with itropo
    ptropo = pclr.isel(layer=itropo.load())

    # mask data and set fill_value pixels to nan
    ptropo = xr.where(tropo_bool, np.nan, ptropo).rename('tropopause_pressure')

    return ptropo


def cal_dphi(saa, vaa):
    '''Calculate the relative azimuth angle (dphi)

    The relative azimuth angle is defined as
        the absolute difference between the viewing azimuth angle and
        the solar azimuth angle.
    It ranges between 0 and 180 degrees
    '''

    dphi = abs(vaa - saa)
    dphi = xr.where(dphi <= 180, dphi, 360 - dphi)

    return dphi.rename('dphi')


def integPr(no2, s5p_p, psfc, ptropo):
    '''
    Integrate the vertical mixing ratios to get vertical column densities

    Because the first TM5 pressure level is the surface pressure
        and we use the tropopause pressure (tp) from TM5,
    we just need to add up the layers from the surface (l = 1)
        up to and including the tropopause level (l = l_{tp}^{TM5}).

    Args:
        no2: no2 vertical mixing ratio (ppp)
        s5p_p: pressure levels (hPa)
        psfc: surface pressure (hPa)
        ptropo: tropopause pressure (hPa)

    Return:
        integrate from psfc to ptropo

    Diagram:
        (Values of "layer" start from 0)

        [tm5_pressure]                           [layer]

        P34(index=33) ++++++++++++++++++++++++++ <--- TM5 top

                      -------------------------- layer32

                      ..........................
                      ..........................

        P26(index=25) ++++++++++++++++++++++++++ <-- Ptropo/integration top

                      -------------------------- layer24 (tropopause)

        P25           ++++++++++++++++++++++++++

                      ..........................
                      ..........................

        P3 (P2 upper) ++++++++++++++++++++++++++

                      -------------------------- layer1 (no2[1])

        P2 (P1 upper) ++++++++++++++++++++++++++

                      -------------------------- layer0 (no2[0])

        P1(index=0)   ++++++++++++++++++++++++++


        For this case, subcolumn is summed from layer0 to layer 32 and
                       sub_layer is from P1 to P26.
        As a result,   vcd should be summed layer0 to layer24,
                       which means sub_layer should be cropped using [:-1, ...]
    '''

    logging.debug(' '*8 + 'Integrating from surface to tropopause ...')

    # constants
    R = 287.3
    T0 = 273.15
    g0 = 9.80665
    p0 = 1.01325e5

    subcolumn = 10 * R * T0 / (g0*p0) \
                   * no2*1e6 \
                   * abs(s5p_p.diff(s5p_p.dims[0])) * 2.6867e16  # DU to moleclues/cm2

    sub_layer = (s5p_p <= psfc) & (s5p_p > ptropo)

    # sum from surface (ground or cloud pressure) to tropopause
    layername = no2.dims[0]
    vcd = subcolumn.where(sub_layer[:-1, ...].values).sum(layername, skipna=False)

    logging.debug(' '*12 + 'Finish integration')

    return vcd


def assign_attrs(da, df_attrs):
    '''assign attributes to DataArray'''
    attrs = df_attrs.loc[da.name]

    return da.assign_attrs(units=attrs.loc['units'],
                           description=attrs.loc['description']
                           )
