# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 11:48:45 2022

@author: Ross
"""
# first of all install h5py library from command line with pip install h5py
# enable h5py
import h5py
import numpy as np

import builtins

#########################################################################################
#################   Change the parameter below  #########################################
#########################################################################################
infile='E:/SpyderProjects3/Prisma/dataset/PRS_L2D_STD_20211103055234_20211103055239_0001/PRS_L2D_STD_20211103055234_20211103055239_0001.he5'
Full=1 # Set value 1 for saving VNIR+SWIR together
if Full:
    outfile='Ahmedabad_full.bin'
    join_priority='SWIR'

else:
    outfile='Ahmedabad_vnir.bin'
    outfile='Ahmedabad_swir.bin'
    
#########################################################################################
#########################################################################################
#########################################################################################


## For image stretching...
from numpy import nanmin, nanmax, where, isnan, fabs
from numpy import nanmean, nanstd

def stddev_stretch(b):
    bf = b.flatten()
    m = nanmean(bf)
    s = nanstd(bf)
    min_ = m - 2*s
    max_ = m + 2*s
    min_ = max(nanmin(bf), min_)
    max_ = min(nanmax(bf), max_)
    return (255.99*(b-min_)/(max_-min_)).clip(0, 255).astype('u1')


def is_string(s):
    return isinstance(s,(str, bytes))

def _write_header_param(fout, paramName, paramVal):
    if paramName.lower() == 'description':
        valStr = '{\n%s}' % '\n'.join(['  ' + line for line
                                       in paramVal.split('\n')])
    elif not is_string(paramVal) and hasattr(paramVal, '__len__'):
        valStr = '{ %s }' % (
            ' , '.join([str(v).replace(',', '-') for v in paramVal]),)
    else:
        valStr = str(paramVal)
    fout.write('%s = %s\n' % (paramName, valStr))


def write_envi_header(fileName, header_dict, is_library=False):
    fout = builtins.open(fileName, 'a+')
    d = {}
    d.update(header_dict)
    # if is_library:
    #     d['file type'] = 'ENVI Spectral Library'
    # elif 'file type' not in d:
    #     d['file type'] = 'ENVI Standard'
    # fout.write('ENVI\n')
    # Write the standard parameters at the top of the file
    std_params = ['description', 'samples', 'lines', 'bands', 'header offset',
                  'file type', 'data type', 'interleave', 'sensor type',
                  'byte order', 'reflectance scale factor', 'map info','bbl']
    for k in std_params:
        if k in d:
            _write_header_param(fout, k, d[k])
    for k in d:
        if k not in std_params:
            _write_header_param(fout, k, d[k])
    fout.close()
# open the PRISMA file
print('Reading Files')
f = h5py.File(infile, 'r')

# reading name and value for root attributes (metadata contained in HDF5 root)
# print('Names of the attributes')
# for attribute in f.attrs:
#     print(attribute,f.attrs[attribute])

# reading names for all attributes (metadata) contained in HDF5 Groups
# specific method for reading the values shall be built depending by the
# specific metadata type (a single value, an array, a matrix, etc)
# def printname(name):
#     print(name)

# f.visit(printname)

# reading SWIR & VNIR datacubes; adjust the “PRS_L1_HCO” string portion
#depending by the specific PRISMA product type (e.g. for L2D product use PRS_L2D_HCO)
    
print('Loading SWIR and VNIR data cubes')
swir = f['/HDFEOS/SWATHS/PRS_L2D_HCO/Data Fields/SWIR_Cube']
vnir = f['/HDFEOS/SWATHS/PRS_L2D_HCO/Data Fields/VNIR_Cube']

# Set this parameter if you want to save data separately or full bands


# list the structure of SWIR data
print("Shape of the data")
print(swir.shape)
# list the structure of VNIR data
print(vnir.shape)

vnir_data=np.zeros(vnir.shape)
swir_data=np.zeros(swir.shape)

### The order is in reverse so we have to sort the order before saving the file

print("Loading Metadata")
wv_vnir=f.attrs['List_Cw_Vnir']
wv_swir=f.attrs['List_Cw_Swir']
wv_vnir=wv_vnir[::-1]
wv_swir=wv_swir[::-1]


fwhm_vnir=f.attrs['List_Fwhm_Vnir']
fwhm_swir=f.attrs['List_Fwhm_Swir']
fwhm_vnir=fwhm_vnir[::-1]
fwhm_swir=fwhm_swir[::-1]


swir.read_direct(swir_data)
vnir.read_direct(vnir_data)


### Transpose the data
# vnir_data=np.transpose(vnir_data,(0,2,1))

vnir_data=np.swapaxes(vnir_data,1,2)
vnir_data=vnir_data[:,:,::-1]

swir_data=np.swapaxes(swir_data,1,2)
swir_data=swir_data[:,:,::-1]


### Scale data for reflectance units
vnir_data=f.attrs['L2ScaleVnirMin']+vnir_data*(f.attrs['L2ScaleVnirMax']-f.attrs['L2ScaleVnirMin'])/65535.0
swir_data=f.attrs['L2ScaleSwirMin']+swir_data*(f.attrs['L2ScaleSwirMax']-f.attrs['L2ScaleSwirMin'])/65535.0


# vnir_data=vnir_data[:,::-1,:]


proj_code=f.attrs['Projection_Id']
proj_name=f.attrs['Projection_Name']
proj_epsg=f.attrs['Epsg_Code']
xmin=min(f.attrs['Product_ULcorner_easting'],f.attrs['Product_LLcorner_easting'])
xmax=max(f.attrs['Product_LRcorner_easting'],f.attrs['Product_URcorner_easting'])
ymin=min(f.attrs['Product_LLcorner_northing'],f.attrs['Product_LRcorner_northing'])
ymax=max(f.attrs['Product_ULcorner_northing'],f.attrs['Product_URcorner_northing'])


from osgeo import ogr, osr

sr=osr.SpatialReference()
sr.ImportFromEPSG(int(proj_epsg))


import matplotlib.pyplot as plt
plt.imshow(stddev_stretch(vnir_data[:,:,[31,19,6]]))

# import data
# data=data.data()

import os
def find_header_file(fname):
    result = fname + '.hdr'
    if os.path.isfile(result):
        return result

    result = os.path.splitext(fname)[0] + '.hdr'
    if os.path.isfile(result):
        return result

    result = fname
    if os.path.isfile(result):
        return result

    return None


from osgeo import gdal

meta={}
meta['wavelength units']='Nanometers'
meta['sensor type']='PRISMA'
meta['data ignore value']= "-9.99e+002"
if not Full:
    
    print('Saving VNIR data')
    
    # Write VNIR File
    rows,cols,bands=vnir_data.shape
    meta['wavelength']=wv_vnir
    meta['fwhm']=fwhm_vnir
    meta['bbl']=[1]*bands
    driver = gdal.GetDriverByName('ENVI')
    outDataset = driver.Create(outfile,
                    cols,rows,bands,gdal.GDT_Float32)
    outDataset.SetProjection(sr.ExportToWkt())
    outDataset.SetGeoTransform((xmin, 30, 0,ymax , 0, -30)) 
    for i in range(bands):
        out = outDataset.GetRasterBand(i + 1)
        out.WriteArray(vnir_data[:, :, i])
        out.FlushCache()
    outDataset = None
    
    hdr_file=find_header_file(outfile)
    write_envi_header(hdr_file, meta)
    
    
    print('Saving SWIR data')
    
    # Write SWIR File
    rows,cols,bands=swir_data.shape
    meta['wavelength']=wv_swir
    meta['fwhm']=fwhm_swir
    meta['bbl']=[1]*bands
    
    driver = gdal.GetDriverByName('ENVI')
    outDataset = driver.Create(outfile,
                    cols,rows,bands,gdal.GDT_Float32)
    outDataset.SetProjection(sr.ExportToWkt())
    outDataset.SetGeoTransform((xmin, 30, 0,ymax , 0, -30)) 
    for i in range(bands):
        out = outDataset.GetRasterBand(i + 1)
        out.WriteArray(swir_data[:, :, i])
        out.FlushCache()
    outDataset = None
    
    hdr_file=find_header_file(outfile)
    write_envi_header(hdr_file, meta)

else:
    
    
    if join_priority=='SWIR':
        full_bands=np.dstack((vnir_data[:,:,np.where(wv_vnir<min(wv_swir))[0]],swir_data ))
        full_wv=list(wv_vnir[np.where(wv_vnir<min(wv_swir))[0]]) + list(wv_swir)
        full_fwhm=list(fwhm_vnir[np.where(wv_vnir<min(wv_swir))[0]]) + list(fwhm_swir)
    else:
        full_bands=np.dstack(vnir_data,swir_data[:,:,np.where(wv_swir>max(wv_vnir))[0]])
        full_wv=list(wv_vnir) + list(wv_swir[np.where(wv_swir>max(wv_vnir))[0]])
        full_fwhm=list(fwhm_vnir)+ list(fwhm_swir[np.where(wv_swir>max(wv_vnir))[0]])
    

    
    print('Saving Full VNIR + SWIR data')
    # Write VNIR File
    rows,cols,bands=full_bands.shape
    
    meta['wavelength']=full_wv
    meta['fwhm']=full_fwhm
    meta['bbl']=[1]*bands
    
    driver = gdal.GetDriverByName('ENVI')
    outDataset = driver.Create(outfile,
                    cols,rows,bands,gdal.GDT_Float32)
    outDataset.SetProjection(sr.ExportToWkt())
    outDataset.SetGeoTransform((xmin, 30, 0,ymax , 0, -30)) 
    for i in range(bands):
        out = outDataset.GetRasterBand(i + 1)
        out.WriteArray(full_bands[:, :, i])
        out.FlushCache()
    outDataset = None
    
    hdr_file=find_header_file(outfile)
    write_envi_header(hdr_file, meta)
