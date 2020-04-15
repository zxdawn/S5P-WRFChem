'''
Some functions for plotting data related to s5p

UPDATE:
    Xin Zhang:
       03/31/2020: Basic

'''

import os
import logging
# import matplotlib
import numpy as np
import pandas as pd
# import seaborn as sns
import proplot as plot
from wrf import ll_to_xy
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def plot_weights_sparse(regridder, output_fig_dir):
    '''Quicklook of the sparse'''
    plt.spy(regridder.weights)
    plt.xlabel('input grid indices')
    plt.ylabel('output grid indices')
    output_name = os.path.join(output_fig_dir, 'weights_sparse.png')
    logging.info(f'Saving quicklook of weights to {output_name}')
    plt.savefig(output_name)

def plot_weights_heatmap(wrf, s5p, regridder, wrf_file, output_fig_dir):
    '''Plot the heatmap of weights for one pixel'''
    # get the value of one pixel around the simulation center
    vcd = s5p['nitrogendioxide_tropospheric_column']
    # assign coords
    vcd = vcd.assign_coords(y=np.arange(vcd.shape[0]),
                            x=np.arange(vcd.shape[1]))

    # select one pixel
    sel_pixel = vcd.where((vcd.longitude >
                           wrf.coords['lon'].mean()-0.1) &
                          (vcd.longitude <
                           wrf.coords['lon'].mean()+0.1) &
                          (vcd.latitude >
                           wrf.coords['lat'].mean()-0.1) &
                          (vcd.latitude <
                           wrf.coords['lat'].mean()+0.1),
                          drop=True)[-1:, -1:]

    # calculate the index for raveled s5p data
    #   because the weight is based on raveled arrays
    pixel_y = sel_pixel.coords['y'] - vcd.coords['y'][0]
    pixel_x = sel_pixel.coords['x'] - vcd.coords['x'][0]
    pixel_indice = pixel_y * len(vcd.coords['x']) + pixel_x

    # get the indices
    match_indice = np.where(np.isin(regridder.weights.row, pixel_indice))
    chem_indice = regridder.weights.col[match_indice]
    # get the x/y in chem grids
    chem_y, chem_x = np.unravel_index(chem_indice,
                                      (wrf['no2'].shape[1], wrf['no2'].shape[2])
                                      )

    df = pd.DataFrame(index=np.unique(chem_y), columns=np.unique(chem_x))
    df = df.fillna(0)

    # fill the df with weights
    for df_index in np.arange(len(chem_y)):
        df.loc[chem_y[df_index], chem_x[df_index]] = \
                regridder.weights.data[match_indice][df_index]

    # proplot version
    f, ax = plot.subplots()
    ax.heatmap(df.columns, df.index,
               df*100, cmap='Blues',
               vmin=0, vmax=4, N=100,
               lw=0.5, edgecolor=None,
               labels=True,
               clip_on=False,
               colorbar='b',
               colorbar_kw={'values': plot.arange(0, 4, 0.5),
                            'label': 'Weights (%)'},
               # labels_kw={'weight': 'bold'},
               )

    # set axis format
    ax.format(alpha=0,
              linewidth=0,
              xticks=[],
              yticks=[],
              # yreverse=True,
              # xloc='top',
              # yloc='right',
              # ticklabelweight='bold',
              # ytickmajorpad=4,
              )

    # convert lon_b/lat_b of the pixel to x/y in chem grid
    nc_ds = Dataset(wrf_file)
    pixel_lon_b = s5p['assembled_lon_bounds'][pixel_y.values[0]:pixel_y.values[0]+2,
                                              pixel_x.values[0]:pixel_x.values[0]+2]
    pixel_lat_b = s5p['assembled_lat_bounds'][pixel_y.values[0]:pixel_y.values[0]+2,
                                              pixel_x.values[0]:pixel_x.values[0]+2]
    x_y = ll_to_xy(nc_ds, pixel_lat_b, pixel_lon_b, as_int=False).values
    # pair together
    x_y = np.stack((x_y[0], x_y[1]), axis=-1)
    x_y[[-2, -1]] = x_y[[-1, -2]]

    # # seaborn version
    # sns.set()
    # f, ax = plt.subplots(figsize=(10, 8))
    # sns.heatmap(df*100,
    #             ax=ax,
    #             cmap="Blues",
    #             # annot_kws={"size": 15},
    #             # cbar_kws={'label': 'weights (%)'},
    #             annot=True, cbar=False,
    #             xticklabels=False,
    #             yticklabels=False)
    # # add % text
    # for t in ax.texts:
    #     t.set_text(t.get_text() + " %")
    # # bug of matplotlib_3.1.1:
    # # https://stackoverflow.com/questions/56942670/
    # #           matplotlib-seaborn-first-and-last-row-cut-in-half-of-heatmap-plot
    # if matplotlib.__version__ == '3.1.1':
    #     bottom, top = ax.get_ylim()
    #     ax.set_ylim(bottom + 0.5, top - 0.5)
    # pair together
    #   we need to convert x/y (chem) to x/y relative to axis,
    #   because the label is at the center, we need to modify values by 0.5
    # x_y = np.stack((x_y[0]-chem_x.min()+0.5, x_y[1]-chem_y.min()+0.5), axis=-1)
    # x_y[[-2, -1]] = x_y[[-1, -2]]

    # add polygon
    poly = Polygon(x_y,
                   edgecolor='orange7',
                   # linewidth=2,
                   fill=None)
    ax.add_patch(poly)

    # save figure
    output_name = os.path.join(output_fig_dir, 'weights_heatmap.png')
    logging.info(f'Saving weights_heatmap to {output_name}')
    f.savefig(output_name)

def plot_regrid_p(wrf, interp_ds, output_fig_dir):
    '''Compare regridded pressure with original WRF pressure'''
    # get variables
    wrf_p = (wrf['P']+wrf['PB'])/100
    regrid_p = interp_ds['p_wrf']

    # subset to model domain
    lon_min = wrf_p.lon.min()
    lon_max = wrf_p.lon.max()
    lat_min = wrf_p.lat.min()
    lat_max = wrf_p.lat.max()
    regrid_p = regrid_p.where((regrid_p.lon >= lon_min) &
                              (regrid_p.lon <= lon_max) &
                              (regrid_p.lat >= lat_min) &
                              (regrid_p.lat <= lat_max),
                              drop=True)

    # set axis
    f, axs = plot.subplots(nrows=1, ncols=2, share=0)
    axs.format(abc=True, abcloc='l', abcstyle='a)')

    # plot surface pressure
    vmin = 950
    vmax = 1000
    wrf_p[0, ...].plot(ax=axs[0],
                       vmin=vmin,
                       vmax=vmax,
                       cbar_kwargs={'label': '(hPa)'})
    regrid_p[0, ...].plot(ax=axs[1],
                          vmin=vmin,
                          vmax=vmax,
                          cbar_kwargs={'label': '(hPa)'})

    # set title
    axs[0].format(title='WRF pressure (hPa)')
    axs[1].format(title='Regridded pressure (hPa)')

    # save figure
    output_name = os.path.join(output_fig_dir, 'comp_regrid.png')
    logging.info(f'Saving comp_regrid to {output_name}')
    f.savefig(output_name)

def plot_interp_p(interp_ds, output_fig_dir):
    '''Compare the interpolated profiles with original simulated profiles'''
    # get the nonan indexes
    # y = 191  # 190
    # x = 312  # 316

    subset = np.where(interp_ds['no2'].notnull())
    # use the first one as the example
    y = subset[1][0]
    x = subset[2][0]

    # set axis
    f, axs = plot.subplots(ncols=2, spanx=0)
    axs.format(abc=True, abcloc='l', abcstyle='a)')

    # no2
    h1 = axs[0].plot(interp_ds['no2'][:, y, x]*1e9,
                     interp_ds['p'][:, y, x],
                     label='TM5')
    h2 = axs[0].plot(interp_ds['no2_wrf'][:, y, x]*1e9,
                     interp_ds['p_wrf'][:, y, x],
                     label='WRF')

    axs[0].legend([h1, h2], loc='ur', frame=False)
    axs[0].format(ylim=(1010, 50),
                  xlim=(0, interp_ds['no2_wrf'][:, y, x].max()*1e9),
                  xlabel='NO$_2$ (ppb)',
                  ylabel='Pressure (hPa)')

    # Temperature
    h1 = axs[1].plot(interp_ds['tk'][:, y, x] - 273.15,
                     interp_ds['p'][:, y, x],
                     label='TM5')
    h2 = axs[1].plot(interp_ds['tk_wrf'][:, y, x] - 273.15,
                     interp_ds['p_wrf'][:, y, x],
                     label='WRF')

    axs[1].legend([h1, h2], loc='ur', frame=False)
    axs[1].format(ylim=(1010, 50),
                  xlim=(-80, 40),
                  xlabel='Temperature ($^\circ$C)',
                  ylabel='Pressure (hPa)')

    # save figure
    output_name = os.path.join(output_fig_dir, 'interp_profile.png')
    logging.info(f'Saving interp_profile to {output_name}')
    f.savefig(output_name)

def plot_bamf_p(bAmfClr, bAmfCld, output_fig_dir):
    for p_index in range(bAmfClr.shape[0]):
        f, axs = plot.subplots(ncols=2)
        x, y = np.meshgrid(bAmfClr.x, bAmfClr.y)
        axs[0].pcolormesh(y, x, bAmfClr[p_index, ...])
        axs[1].pcolormesh(y, x, bAmfCld[p_index, ...])

        # save figure
        output_name = os.path.join(output_fig_dir, f'bAmf_layer{p_index}.png')
        logging.info(f'Saving bAmf to {output_name}')
        f.savefig(output_name)
