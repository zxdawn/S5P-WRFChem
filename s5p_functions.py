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
        """ Reads in the settings.txt file

        Example 'settings.txt':
        [ALL]
        start_date = 2019-07-25     ; The beginning of s5p data (yyyy-mm-dd)
        end_date = 2019-07-25       ; The end of s5p data (yyyy-mm-dd)
        s5p_nc_dir = s5p_data       ; where s5p data are saved
        wrf_nc_dir = wrf_data       ; where wrfchem data are saved
        output_data_dir = output    ; where generated files are saved
        output_fig_dir = figures    ; where figures are saved
        overwrite = False           ; whether overwrite the same generated files (True/False)
        plot_weights = False        ; whether plot the weights
        """

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
                          input_p_dimname: input_array.shape[0],
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

    # interpolate data and keep the original order of dims
    interp_no2 = xr_interp(regrid_vars['no2'],
                           regrid_vars['p'], 'bottom_top',
                           s5p['p'].load(), 'layer').transpose(
                           'layer',
                           ...,
                           transpose_coords=False)

    interp_tk = xr_interp(regrid_vars['tk'],
                          regrid_vars['p'], 'bottom_top',
                          s5p['p'].load(), 'layer').transpose(
                          'layer',
                          ...,
                          transpose_coords=False)

    # save all variables into one Dataset
    interp_ds = xr.Dataset({
                            'no2': interp_no2,
                            'tk': interp_tk,
                            'p': s5p['p'],
                            'no2_wrf': regrid_vars['no2'],
                            'tk_wrf': regrid_vars['tk'],
                            'p_wrf': regrid_vars['p'],
                           })

    logging.info(' '*8 + 'Finish interpolation')

    return interp_ds

def interp_to_cld(profile, pclr, pcld):
    '''Interpolate data to pressure levels including cloud pressure'''
    return xr_interp(profile,
                     pclr, 'layer',
                     pcld, 'plevel').transpose(
                     'plevel',
                     ...,
                     transpose_coords=False)

def concat_p(pcld, ptm5):
    '''Concatenate cloud pressures and tm5 pressures'''
    # Check submitted issue by Xin:
    #   https://github.com/pydata/xarray/issues/3954
    pcld = pcld.expand_dims(layer=[ptm5.sizes['layer']]).drop(['layer'])
    s5p_pcld = xr.concat([ptm5, pcld], dim='layer')

    # sort pressures in DataFrame
    # sort_values isn't available for dask DataFrame yet
    #   https://github.com/dask/dask/issues/958
    logging.info(' '*6 + 'Sorting pressures slowly...')
    df = s5p_pcld.to_series()
    df = df.reset_index() \
           .sort_values(['y', 'x', 'tm5_pressure'], ascending=False) \
           .set_index(['layer', 'y', 'x'])

    # Unfortunately, `set_levels` has bug ...
    # Check issues submitted by Xin:
    #   https://github.com/pandas-dev/pandas/issues/33420
    #   https://github.com/pydata/xarray/issues/3957
    # df.index.set_levels(np.arange(s5p_pcld.sizes['layer']),
    #                     level='layer',
    #                     inplace=True)

    # set new indexes for 'layer'
    df.index = pd.MultiIndex.from_arrays(
                    [np.tile(np.arange(s5p_pcld.sizes['layer']),
                             s5p_pcld.sizes['y']*s5p_pcld.sizes['x']),
                     df.index.get_level_values(1),
                     df.index.get_level_values(2)]) \
                 .set_names('plevel', level=0)

    # convert DataFrame to dask DataArray
    s5p_pcld = df.to_xarray()['tm5_pressure'] \
                 # .chunk({'plevel': s5p_pcld.shape[0],
                 #         'y': s5p_pcld.shape[1],
                 #         'x': s5p_pcld.shape[2]}) \

    # we need to drop coords for y and x first
    #   https://github.com/pydata/xarray/issues/2283
    s5p_pcld = s5p_pcld.drop(['y', 'x'])

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
    # mask with where and set fill_value pixels to nan
    ptropo = xr.where(tropo_bool, np.nan, ptropo).rename('tropopause_pressure')

    return ptropo

def cal_rza(sza, vza):
    '''Calculate the relative azimuth angle

    The relative azimuth angle is defined as
        the absolute difference between the satellite azimuth angle and
        the solar azimuth angle.
    It ranges between 0 and 180 degrees

    '''

    rza = abs(sza - vza)
    rza = xr.where(rza > 180, 360 - rza, rza).rename('dphi')

    return rza

def integPr(ds):
    '''
    Integrate the vertical mixing ratios to get vertical column densities

    Because the first TM5 pressure level is the surface pressure
        and we use the tropopause pressure (tp) from TM5,
    we just need to add up the layers from the surface (l = 1)
        up to and including the tropopause level (l = l_{tp}^{TM5}).

    Args:
        vmr: vertical mixing ratio (ppp)
        p: pressure levels (hPa)
        psfc: surface pressure (hPa)
        ptropo: tropopause pressure (hPa)
        layername: the name of pressure layer

    Return:
        integrate up to tropopause pressure

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

        P3            ++++++++++++++++++++++++++

                      -------------------------- layer1

        P2            ++++++++++++++++++++++++++

                      -------------------------- layer0

        P1(index=0)   ++++++++++++++++++++++++++


        For this case, subcolumn is summed from layer0 to layer 32 and
                       sub_layer is from P1 to P26.
        As a result,   vcd should be summed layer0 to layer24,
                       which means sub_layer should be cropped using [1:, ...]
    '''

    logging.debug(' '*8 + 'Integrating from surface to tropopause ...')

    # molecular mass (kg)     * g (m/s2) * (Pa/hPa) * (m2/cm2)
    mg = (28.97/6.02e23)*1e-3 * 9.8 * 1e-2 * 1e4
    # subcolumn[i, :, :] = (vmr[i]+vmr[i+1]) * (p[i]-p[i+1]) / (2*mg)

    layername = ds['no2'].dims[0]
    subcolumn = ds['no2'].rolling({layername: 2}).sum()[1:, ...] \
        * abs(ds['tm5_pressure'].diff(ds['tm5_pressure'].dims[0])) / (2*mg)

    # get the layer mask
    sub_layer = (ds['tm5_pressure'] <= ds['surface_pressure']) & \
        (ds['tm5_pressure'] >= ds['tropopause_pressure'])

    # sum from surface (ground or cloud pressure) to tropopause
    vcd = subcolumn.where(sub_layer[1:, ...]).sum(layername)

    # set 0 to nan
    vcd = vcd.where(vcd > 0)

    logging.debug(' '*12 + 'Finish integration')

    return vcd

def assign_attrs(da, df_attrs):
    attrs = df_attrs.loc[da.name]
    return da.assign_attrs(units=attrs.loc['units'],
                           description=attrs.loc['description']
                           )
