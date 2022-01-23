'''
Some utils for s5p_main.py

UPDATE:
    Xin Zhang:
       03/31/2020: Basic
       04/01/2020: Add the interp part
       04/12/2020: Finish the AMF part
       03/09/2021: Fix the AMF bug
'''

import os
import glob
import ntpath
import logging
import numpy as np
import xesmf as xe
import xarray as xr
import pandas as pd
from xgcm import Grid
from satpy.scene import Scene
from datetime import datetime
from xgcm.autogenerate import generate_grid_ds
from s5p_functions import *


def validation(cfg):
    '''Validate inputs'''
    validate_date(cfg['start_date'])
    validate_date(cfg['end_date'])
    validate_path(cfg['s5p_nc_dir'], 's5p_nc_dir')
    validate_path(cfg['wrf_nc_dir'], 'wrf_nc_dir')
    validate_path(cfg['output_data_dir'], 'output_data_dir')
    validate_path(cfg['output_fig_dir'], 'output_fig_dir')


def load_s5p(date_in, s5p_nc_dir, tm5_prof):
    '''Load s5p data'''
    # get the s5p data by the datetime
    f_s5p_pattern = os.path.join(s5p_nc_dir, f'*{date_in.strftime("%Y%m%d")}*')
    f_s5p = glob.glob(f_s5p_pattern)
    logging.info(' '*4+f'Reading {f_s5p} ...')
    s5p = Scene(f_s5p, reader='tropomi_l2')

    vnames = ['nitrogendioxide_slant_column_density',
              'nitrogendioxide_stratospheric_column',
              'nitrogendioxide_tropospheric_column',
              'nitrogendioxide_ghost_column',
              'assembled_lat_bounds', 'assembled_lon_bounds',
              'latitude', 'longitude',
              'latitude_bounds', 'longitude_bounds',
              'surface_albedo_nitrogendioxide_window', 'surface_pressure',
              'cloud_pressure_crb', 'cloud_albedo_crb',
              'cloud_fraction_crb_nitrogendioxide_window',
              'cloud_radiance_fraction_nitrogendioxide_window',
              'solar_azimuth_angle', 'viewing_azimuth_angle',
              'solar_zenith_angle', 'viewing_zenith_angle',
              'tm5_constant_a', 'tm5_constant_b',
              'tm5_tropopause_layer_index',
              'qa_value',
              'no2_scd_flag',  # only available for processed cloudy data
              'time_utc', 'air_mass_factor_troposphere',
              'air_mass_factor_clear', 'air_mass_factor_cloudy', 'amf_geo',
              'air_mass_factor_stratosphere']

    if tm5_prof:
        # Read the TM5-MP a-priori profile
        vnames.extend(['no2_vmr', 'temperature'])

    logging.debug(' '*4 + f'Reading vnames: {vnames}')
    s5p.load(vnames)
    # another option: load all available variables
    # s5p.load(s5p.all_dataset_names())

    # using pandas to convert string to timestamp, then to datetime without tz.
    mean_t = pd.to_datetime(s5p['time_utc'].values).mean() \
        .to_pydatetime().replace(tzinfo=None)

    # set global attrs
    s5p.attrs['s5p_filename'] = os.path.basename(f_s5p[0])

    # calculate pressure levels
    a = s5p['tm5_constant_a']
    b = s5p['tm5_constant_b']
    psfc = s5p['surface_pressure']

    low_p = (a[:, 0] + b[:, 0]*psfc)/1e2
    high_p = (a[:, 1] + b[:, 1]*psfc)/1e2

    s5p['p'] = xr.concat([low_p, high_p.isel(layer=-1)], dim='layer')
    s5p['p'] = s5p['p'].rename('tm5_pressure')
    s5p['p'].attrs['units'] = 'hPa'

    # read lut
    lut_pattern = os.path.join(s5p_nc_dir, 'S5P_OPER_LUT_NO2AMF*')
    lut = xr.open_mfdataset(lut_pattern, combine='by_coords')

    logging.info(' '*8 + 'Finish reading')

    return s5p, vnames, lut, mean_t


def load_wrf(date_in, wrf_nc_dir):
    '''Load wrf data

    Because loading large files could cost much time,
    please use "frames_per_outfile=1" in "namelist.input"
    to generate wrfout* files.

    '''

    # get all wrf output files in the same day
    f_wrf_pattern = os.path.join(wrf_nc_dir, date_in.strftime('wrfout_*_%Y-%m-%d_*'))
    wrf_list = glob.glob(f_wrf_pattern)

    # omit the directory and get "yyyy-mm-dd_hh:mm:ss"
    wrf_dates = [datetime.strptime(
                            ntpath.basename(name).split('_', maxsplit=2)[-1],
                            '%Y-%m-%d_%H:%M:%S')
                 for name in wrf_list]

    # get the index of the closest datetime
    wrf_index = wrf_dates.index(min(wrf_dates, key=lambda d: abs(d - date_in)))
    wrf_file = wrf_list[wrf_index]

    # read selected wrf data
    logging.info(' '*4 + f'Reading {wrf_file} ...')
    wrf = xr.open_dataset(wrf_file)
    wrf.attrs['wrfchem_filename'] = os.path.basename(wrf_file)

    # generate lon_bounds and lat_bounds
    wrf = wrf.rename({'XLONG': 'lon', 'XLAT': 'lat'}).isel(Time=0)
    wrf = generate_grid_ds(wrf, {'X': 'west_east', 'Y': 'south_north'},
                           position=('center', 'outer'))

    # convert from potential temperature to absolute temperature (K)
    # http://mailman.ucar.edu/pipermail/wrf-users/2010/001896.html
    # http://mailman.ucar.edu/pipermail/wrf-users/2013/003117.html
    wrf['T'] = (wrf['T'] + 300) * ((wrf['P'] + wrf['PB']) / 1e5) ** 0.2865

    # generate grid and bounds
    grid_ds = Grid(wrf, periodic=False)
    bnd = 'extrapolate'
    wrf.coords['lon_b'] = grid_ds.interp(grid_ds.interp(wrf['lon'], 'X',
                                                        boundary=bnd,
                                                        fill_value=np.nan),
                                         'Y',
                                         boundary=bnd,
                                         fill_value=np.nan)
    wrf.coords['lat_b'] = grid_ds.interp(grid_ds.interp(wrf['lat'], 'X',
                                         boundary=bnd,
                                         fill_value=np.nan),
                                         'Y',
                                         boundary=bnd,
                                         fill_value=np.nan)
    logging.info(' '*8 + 'Finish reading')

    return wrf_file, wrf


def regrid(wrf, s5p, tm5_prof):
    '''Get bounds and regrid data to s5p pixels

    Currently, the boundary for WRF-Chem option is wrong:
        https://github.com/JiaweiZhuang/xESMF/issues/17

    Temporary solution:
        https://github.com/JiaweiZhuang/xESMF/issues/15
        https://github.com/JiaweiZhuang/xESMF/issues/22
    '''

    if tm5_prof:
        logging.info(' '*4 + 'Regridding TM5 data directly ...')
        interp_ds = xr.Dataset({
                                'no2': s5p['no2_vmr'].transpose('layer',
                                                                ...,
                                                                transpose_coords=False),
                                'tk': s5p['temperature'].transpose('layer',
                                                                   ...,
                                                                   transpose_coords=False),
                               })

        regridder = None

    else:
        logging.info(' '*4 + 'Regridding WRF-Chem data to pixels ...')

        # create bounds of wrf and s5p
        wrf_bounds = {'lon': wrf.lon,
                      'lat': wrf.lat,
                      'lon_b': wrf.lon_b,
                      'lat_b': wrf.lat_b,
                      }

        s5p_bounds = {'lon': s5p['longitude'],
                      'lat': s5p['latitude'],
                      'lon_b': s5p['assembled_lon_bounds'],
                      'lat_b': s5p['assembled_lat_bounds'],
                      }

        # create regridder
        regrid_method = 'conservative_normed'
        regridder = xe.Regridder(wrf_bounds, s5p_bounds,
                                 regrid_method,
                                 # filename=f'{regrid_method}_wrf_s5p.nc',
                                 # reuse_weights=True,
                                 )

        # create regridder
        regrid_vars = {'no2': regridder(wrf['no2'])/1e6,  # ppm --> ppp
                       'o3': regridder(wrf['o3'])/1e6,  # ppm --> ppp
                       'tk': regridder(wrf['T']),
                       'p': regridder(wrf['P']+wrf['PB'])/100,  # hPa
                       }

        # apply regridder
        if 'no' in wrf.keys():
            # NO is for the retrieval of LNOx
            # It's fine if you don't care about the LNOx
            regrid_vars.update({'no': regridder(wrf['no'])/1e6})  # ppm --> pp
        if 'o3' in wrf.keys():
            regrid_vars.update({'o3': regridder(wrf['o3'])/1e6})  # ppm --> pp

        # set 0 to nan in regridded variables
        #   as no simulation data is available for these pixels
        for key in regrid_vars:
            regrid_vars[key] = regrid_vars[key].where(regrid_vars[key] > 0)

        logging.info(' '*8 + 'Finish Regridding')

        # save to Dataset
        interp_ds = interp_to_tm5(regrid_vars, s5p)

    return regridder, interp_ds


def cal_bamf(s5p, lut):
    '''Calculate the Box-AMFs based on the LUT file

    Args:
        - albedo: Surface Albedo
        - dphi: Relative azimuth angle
        - mu: Cosine of viewing zenith angle
        - mu0: Cosine of solar zenith angle
        - p: Pressure Levels
        - p_surface: surface_air_pressure
        - sza: Solar zenith angle
        - vza: Viewing zenith angle
        - amf: Box air mass factor

    '''

    logging.info(' '*4 + 'Calculating box-AMFs using LUT ...')

    new_dim = ['y', 'x']

    # get vars from s5p data
    albedo = xr.DataArray(s5p['surface_albedo_nitrogendioxide_window'],
                          dims=new_dim)
    cloud_albedo = xr.DataArray(s5p['cloud_albedo_crb'],
                                dims=new_dim)

    # use surface pressure from TROPOMI (input from ERA-Interim)
    p_surface = xr.DataArray(s5p['surface_pressure']/1e2,  # hPa
                             dims=new_dim)
    p_cloud = xr.DataArray(s5p['cloud_pressure_crb']/1e2,  # hPa
                           dims=new_dim)

    # calculate angles
    dphi = xr.DataArray(cal_dphi(s5p['solar_azimuth_angle'],
                                 s5p['viewing_azimuth_angle']),
                        dims=new_dim)
    mu0 = xr.DataArray(np.cos(np.deg2rad(s5p['solar_zenith_angle'])),
                       dims=new_dim)
    mu = xr.DataArray(np.cos(np.deg2rad(s5p['viewing_zenith_angle'])),
                      dims=new_dim)

    da = lut['amf'].assign_coords(p=np.log(lut['amf'].p), p_surface=np.log(lut['amf'].p_surface))
    # da = da.where(da>0)

    # interpolate data by 2d arrays
    '''
    if you meet "KeyError: nan", see the method below:
        although xarray >= 0.16.1 fix this issue,
        regrid_dataset broken with xarray=0.16.1
            (https://github.com/pangeo-data/xESMF/pull/47)
        so, check https://github.com/pydata/xarray/pull/3924
            and edit the xarray/core/missing.py file by yourself
    '''
    bAmfClr_p = da.interp(albedo=albedo.clip(0,1),
                          p_surface=np.log(p_surface),
                          dphi=dphi,
                          mu0=mu0,
                          mu=mu)

    bAmfCld_p = da.interp(albedo=cloud_albedo.clip(0,1),
                          p_surface=np.log(p_cloud),
                          dphi=dphi,
                          mu0=mu0,
                          mu=mu)

    # interpolate to TM5 pressure levels
    bAmfClr = xr_interp(bAmfClr_p,
                        bAmfClr_p.coords['p'], 'p',
                        np.log(s5p['p'].rolling({'layer': 2}).mean()[1:, ...].load()), 'layer').transpose(
                        'layer',
                        ...,
                        transpose_coords=False)

    bAmfCld = xr_interp(bAmfCld_p,
                        bAmfCld_p.coords['p'], 'p',
                        np.log(s5p['p'].rolling({'layer': 2}).mean()[1:, ...].load()), 'layer').transpose(
                        'layer',
                        ...,
                        transpose_coords=False)

    # because the bAMF is normalized by amf_geo in the LUT file, we need to multiply bAMFs by amf_geo
    bAmfClr *= s5p['amf_geo']
    bAmfCld *= s5p['amf_geo']

    # convert s5p['p'] to dask array
    s5p['p'] = s5p['p'].chunk({'layer': s5p['p'].shape[0],
                               'y': s5p['p'].shape[1],
                               'x': s5p['p'].shape[2]})

    logging.info(' '*8 + 'Finish calculating box-AMFs')

    return bAmfClr, bAmfCld, [albedo, p_surface, cloud_albedo, p_cloud, dphi, mu0, mu]


def cal_amf(s5p, interp_ds, bAmfClr, bAmfCld, del_lut):
    '''Calculate AMFs'''
    logging.info(' '*4 + 'Calculating AMFs based on box-AMFs ...')

    # get simulated profiles
    no2 = interp_ds['no2']
    tk = interp_ds['tk']

    # for LNOx research
    if 'no' in interp_ds.keys():
        no = interp_ds['no']
    if 'o3' in interp_ds.keys():
        o3 = interp_ds['o3']

    # the temperature correction factor, see TROPOMI ATBD file
    ts = 220  # temperature of cross-section [K]
    factor = 1 - 0.00316*(tk-ts) + 3.39e-6*(tk-ts)**2

    # load variables
    psfc = s5p['surface_pressure'] / 1e2  # hPa
    pcld = s5p['cloud_pressure_crb'] / 1e2  # hPa
    cf = s5p['cloud_fraction_crb_nitrogendioxide_window']
    crf = s5p['cloud_radiance_fraction_nitrogendioxide_window']
    itropo = s5p['tm5_tropopause_layer_index']
    s5p_pclr = s5p['p']
    ptropo = cal_tropo(s5p_pclr, itropo)

    # set units
    psfc.attrs['units'] = 'hPa'
    pcld.attrs['units'] = 'hPa'

    # concatenate surface pressures, cloud pressures and tm5 pressures
    s5p_pcld = concat_p([psfc, pcld], s5p_pclr)

    # get the scattering weights
    bAmfClr = bAmfClr * factor
    bAmfCld = bAmfCld * factor

    # interpolate profiles to pressure levels including cloud pressure
    no2 = interp_to_layer(no2, s5p_pclr, s5p_pcld)
    bAmfClr = interp_to_layer(bAmfClr, s5p_pclr, s5p_pcld)
    bAmfCld = interp_to_layer(bAmfCld, s5p_pclr, s5p_pcld)
    clearSW = no2 * bAmfClr
    cloudySW = no2 * bAmfCld

    # for LNOx research
    if 'no' in interp_ds.keys():
        no = interp_to_layer(no, s5p_pclr, s5p_pcld).rename('noapriori')
    if 'o3' in interp_ds.keys():
        o3 = interp_to_layer(o3, s5p_pclr, s5p_pcld).rename('o3apriori')

    # logging.info(' '*6 + 'Calculating ghost column ...')
    # ghost = integPr(no2, s5p_pcld, psfc, pcld.rename('tropopause_pressure')).rename('ghost_column')

    # integrate from surface pressure to tropopause
    logging.info(' '*6 + 'Calculating vcdGnd ...')
    vcdGnd = integPr(no2, s5p_pcld, psfc, ptropo).rename('vcdGnd')

    # integrate from cloud pressure to tropopause
    logging.info(' '*6 + 'Calculating vcdCld ...')
    vcdCld = integPr(no2, s5p_pcld, pcld, ptropo).rename('vcdCld')

    # for LNOx research
    if 'no' in interp_ds.keys():
        logging.info(' '*6 + 'Calculating vcdGnd_no ...')
        vcdGnd_no = integPr(no, s5p_pcld, psfc, ptropo).rename('vcdGnd_no')

    logging.info(' '*6 + 'Calculating scdClr ...')
    scdClr = integPr(clearSW, s5p_pcld, psfc, ptropo).rename('scdClr')

    logging.info(' '*6 + 'Calculating scdCld ...')
    scdCld = integPr(cloudySW, s5p_pcld, pcld, ptropo).rename('scdCld')

    # #set Cld DataArays to nan again for "clear" pixels
    # cld_pixels = (cf > 0) & (crf > 0) & (pcld > ptropo)
    # vcdCld = vcdCld.where(cld_pixels, 0)
    # scdCld = scdCld.where(cld_pixels, 0)

    # calculate AMFs
    logging.info(' '*6 + 'Calculating amfClr ...')
    amfClr = scdClr / vcdGnd
    # https://github.com/pydata/xarray/issues/2283
    # amfClr = amfClr.where((crf != 1) & (cf != 1), 0)

    logging.info(' '*6 + 'Calculating amfCld ...')
    amfCld = scdCld / vcdGnd
    # amfCld = amfCld.where(cld_pixels, 0)

    logging.info(' '*6 + 'Calculating amf and amfVis ...')
    amf = amfCld*crf + amfClr*(1-crf)
    amfVis = amf*vcdGnd / (vcdCld*cf+vcdGnd*(1-cf))

    # calculate averaging kernel
    logging.info(' '*6 + 'Calculating averaging kernel ...')
    sc_weights = crf*bAmfCld.where((s5p_pcld.rolling({s5p_pcld.dims[0]: 2}).mean()[1:, ...] < ptropo), 0) \
        + (1-crf)*bAmfClr
    avgKernel = sc_weights / amf

    # calculate vertical column densities
    scdTrop = s5p['nitrogendioxide_tropospheric_column'] * s5p['air_mass_factor_troposphere']
    no2Trop = scdTrop / amf
    no2TropVis = scdTrop / amfVis

    # rename DataArrays
    amf = amf.rename('amfTrop')
    amfVis = amfVis.rename('amfTropVis')
    bAmfClr = bAmfClr.rename('swClr')
    bAmfCld = bAmfCld.rename('swCld')
    avgKernel = avgKernel.rename('avgKernel')
    no2 = no2.rename('no2apriori')
    s5p_pcld = s5p_pcld.rename('plevels')
    no2Trop = no2Trop.rename('no2Trop')
    no2TropVis = no2Trop.rename('no2TropVis')

    # read attributes table
    df_attrs = pd.read_csv('attrs_table.csv',
                           sep=' *, *',  # delete spaces
                           engine="python"
                           ).set_index('name')

    # drop useless coordinates and assign attributes
    da_list = []
    saved_da = [amf, amfVis, bAmfClr, bAmfCld, avgKernel,
                no2, no2Trop, no2TropVis,
                scdClr, scdCld, vcdGnd,
                ptropo, s5p_pcld.isel(plevel=slice(None, -1))]

    if 'no' in interp_ds.keys():
        if 'o3' in interp_ds.keys():
            saved_da.extend([o3, no, vcdGnd_no])
        else:
            saved_da.extend([no, vcdGnd_no])

    for da in saved_da:
        da_list.append(assign_attrs(da.drop(list(da.coords)), df_attrs))

    # merge to one Dataset
    ds = xr.merge(da_list)

    return ds


def save(s5p, wrf, ds, vnames, output_file, tm5_prof):
    '''Save data to nc files'''
    # it is better to load before save to netcdf
    # https://github.com/pydata/xarray/issues/2912

    # set global attributes
    header_attrs = {'wrfchem_filename': wrf.attrs['wrfchem_filename'],
                    's5p_filename': s5p.attrs['s5p_filename'],
                    }

    # set compression
    comp = dict(zlib=True, complevel=9)

    # save loaded s5p data
    logging.info(' '*4 + f'Saving to {output_file}')

    s5p.save_datasets(filename=output_file,
                      datasets=vnames,
                      groups={'S5P': vnames},
                      compute=True,
                      header_attrs=header_attrs,
                      writer='cf',
                      engine='netcdf4',
                      compression=comp,
                      )

    # save data into different groups
    if tm5_prof:
        group_name = 'TM5'
    else:
        group_name = 'CHEM'

    # save to one netcdf file
    encoding = {var: comp for var in ds.data_vars}
    ds.load().to_netcdf(path=output_file,
                        mode='a',
                        engine='netcdf4',
                        group=group_name,
                        encoding=encoding)
