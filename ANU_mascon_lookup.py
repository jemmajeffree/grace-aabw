import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.spatial
import xarray as xr

def lookup_mascon(mascon_definition_filepath, 
                  lonlat):
        
    '''lonlat is an array with shape (2, ...), although a 1D tuple or list should also work'''

    def to_xyz(lonlat):
        '''Turns longitude and latitude coordinates into x,y,z coordinates
        Spherical earth is sufficient for determining the nearest ternary mascon
        This is just a way of making sure I can wrap around 360 and that nearest longitude 
        doesn't dominate like it would at the poles in euclidian lat/lon distance

        lonlat is an array with shape (2, ...), although a 1D tuple or list should also work'''

        #Using trig identity cheats halves the time this function takes - possibly worth implementing but be careful
        #of Nans from sqrt(-0.0001) type things

        lon = lonlat[0]
        lat = lonlat[1]

        #Attempt to make sure lat & lon are the right way around
        assert np.all(np.abs(lat)<91), 'lat outside of range: '+str(lat[np.abs(lat)>91])
        R_E = 6.371*10**6 #Earth's radius

        z = np.sin(np.deg2rad(lat))  *R_E
        xy = np.cos(np.deg2rad(lat))  *R_E
        x = np.cos(np.deg2rad(lon))*xy #I'm not certain these are conventionally defined (x & y might be swapped), 
        y = np.sin(np.deg2rad(lon))*xy #but it doesn't matter for this application.

        return np.array((x,y,z))
    
    print('Building mascon lookup tree (~60 seconds)')    
    #Build primary lookup for any x,y,z (15 seconds)

    ternary_list = [] # If you know exactly how many ternary mascons there are 
    primary_list = [] # and wanted to make this an array it would go faster

    with open(mascon_definition_filepath) as inputFile:
        while True:
            l = inputFile.readline().split()
            if l ==[]: #Stop at the end of the file
                break
            if l[0] == '#': #Skip commented out lines
                continue

            if l[1][0] == 'T': #I only care about the ternary mascons
                primary = int(l[9])
                lat = round(float(l[2]),5)
                lon = round(float(l[3])%360,5)

                ternary_list.append((lon,lat))
                primary_list.append(primary)   #Indices for these two lists should correlate

    ternary_xyz = to_xyz(np.array(ternary_list).T).T
    ternary_tree = scipy.spatial.KDTree(ternary_xyz)
    
    #Build a dictionary primary lookup of x,y,z points (15 seconds)
    primary_name_lookup = {}
    for i in range(len(ternary_list)):
        key = tuple(np.round(ternary_xyz[i],6))
        primary_name_lookup[key] = primary_list[i]
    
    
    #Build lookup for each element in the tree (20 seconds)
    def what_primary(xyz):
        ''' return what primary mascon the xyz coordinates specified belong to'''
        key = tuple(np.round(xyz.squeeze(),6))
        return primary_name_lookup[key]

    tree_to_primary = np.zeros(ternary_tree.data.shape[0])
    for i,d in enumerate(ternary_tree.data):
        tree_to_primary[i] = what_primary(d)
        
  
    input_xyz = to_xyz(lonlat)
    
    #Actual heavy lifting. Possibly 30 minutes
    print('Starting lookup for supplied indices. This could take a while \n \
(~2 min? It took 30 min once but I haven\'t been able to replicate that)')
    distances, indexes = ternary_tree.query(input_xyz.T)
    print('Finished mascon lookup')
    return tree_to_primary[indexes]

def lookup_mascon_grid(mascon_definition_filepath, 
                  lon_coords, #array-like of longitude (shape M)
                  lat_coords, #array-like of latitude (shape N)
                  output_filepath, #I don't ever want to calculate that huge a thing without dumping to netcdf
                 ):
    
    grid = np.array(np.meshgrid(lon_coords,lat_coords, indexing='ij'))
    flat_grid = grid.reshape(2,lon_coords.shape[0]*lat_coords.shape[0])
    primary_grid = lookup_mascon(mascon_definition_filepath, grid).reshape(lon_coords.shape[0],lat_coords.shape[0])
    
    #Save output
    output = xr.Dataset({'primary_mascon':xr.DataArray(primary_grid,dims=('xt_ocean','yt_ocean'),coords={'xt_ocean':lon_coords,'yt_ocean':lat_coords})})
    output.to_netcdf(output_filepath)
            
    return output
    
    
    
    