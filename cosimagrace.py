import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import sklearn.metrics
import scipy.stats
import matplotlib as mpl

def transect_expt_id(lat,lon1,lon2,expt='',freq='1 monthly',min_rho = 1037):
    '''Make a unique id string for this setup'''
    id_str = str(lat)+'N_'+str(lon1)+'_'+str(lon2)+'E'
    
    if not(freq == '1 monthly'):
        id_str+='_'+freq
    
    if expt == '01deg_jra55v13_ryf9091':
        id_str += '_ryf'
    elif expt =='01deg_jra55v140_iaf_cycle3':
        id_str += '_iaf'
    else:
        id_str+='_'+expt
        
    if not(min_rho == 1037):
        id_str+='_'+str(round(minrho,4))
        
    return id_str

def nc_name(lat, lon1,lon2,expt='',freq='1 monthly',min_rho = 1037):
    '''define where the data I'm looking for is. Copied across from Calculating_overturning.ipynb'''
    name = '/g/data/x77/jj8842/calculated_overturning/overturning_'+transect_expt_id(lat,lon1,lon2,expt,freq,min_rho)+'.nc'     
    return name

def plot_transects(t):
    '''Plot as lines on a map where the transects are
    t = xarray object of transects'''
    ### FYI, I think there's problems with this code if transect.id is not 0,1,2...
    
    
    for i in t.id:
        l = np.zeros((2,2))
        l[0] = t.sel(id=i).lon
        l[1] = t.sel(id=i).lat
        
        if np.isnan(l[0,0]):
            l[0,0] = -285
        
        if np.isnan(l[0,1]):
            l[0,1] = 85

        if l[0,1]>l[0,0]:
            plt.plot(l[0],l[1],c='k',marker='+',markersize=10,linewidth=1)
        else:
            plt.plot((-285,l[0,1]),l[1],c='k',marker='+',markersize=10,linewidth=1)
            plt.plot((l[0,0],85),l[1],c='k',marker='+',markersize=10,linewidth=1)
            
def plot_transects_fancy(t,color,width,cmap,vmin = None, vmax = None, outline_t = [],extend=None):
    ### FYI, I think there's problems with this code if transect.id is not 0,1,2...
    '''Plot as lines on a map where the transects are
    t = xarray object of transects'''
    
    if vmin is None:
        vmin = np.min(color)
    if vmax is None:
        vmax = np.max(color)
    for i,id in enumerate(t.id):
        c = cmap((color[i]-vmin)/(vmax-vmin))
        if id in outline_t:
            edgecolor='k'
        else:
            edgecolor='None'
        l = np.zeros((2,2))
        l[0] = t.sel(id=id).lon
        l[1] = t.sel(id=id).lat

        if np.isnan(l[0,0]):
            l[0,0] = -285

        if np.isnan(l[0,1]):
            l[0,1] = 85

        if l[0,1]>l[0,0]:
            plt.fill_between(l[0],l[1]+width/2,l[1]-width/2,facecolor=c,edgecolor=edgecolor)
        else:
            plt.fill_between((-285,l[0,1]),l[1]+width/2,l[1]-width/2,facecolor=c,edgecolor=edgecolor)
            plt.fill_between((l[0,0],85),l[1]+width/2,l[1]-width/2,facecolor=c,edgecolor=edgecolor)
            
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    return plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),extend=extend)
            
def plot_over_transects(R2, transect, stringify = lambda x : str(np.round(x,2))):
    ### FYI, I think there's problems with this code if transect.id is not 0,1,2...
    '''Plot the numbers R2 belonging to the transects ti
    stringify turns floats (or whatever is in R2) into strings for plotting'''
    for i,e in enumerate(transect.id):
        lat = np.mean(transect.lat[e])-4
        lon = np.mean(transect.lon[e])
        if np.isnan(lon):
            lon = -250
        elif transect.lon[e][0]>transect.lon[e][1]:
            lon = transect.lon[e][1]-10
        plt.text(lon,lat,stringify(R2[i]),c='k',fontsize=16,horizontalalignment='center')
        
def sel_longitude(data,lon1,lon2,coord='grid_xt_ocean'):
    '''Select a longitude range of data. If this longitude range loops around the -280 80 gap, roll data before selecting'''
    
    if np.isnan(lon1):
        lon1 = -280
    if np.isnan(lon2):
        lon2 = 80
    if lon1<lon2:
        return data.sel({coord:slice(lon1,lon2)})
    
    n_roll = int(np.sum(data.coords[coord] >lon1))
    newdata    =  data.roll({coord:n_roll},roll_coords=True)
    
    new_coords =  np.array(newdata.coords[coord])
    new_coords[np.where(new_coords>float(lon1))]-=360
    newdata = newdata.assign_coords({coord:new_coords})
    
    return newdata.sel({coord:slice(lon1-360,lon2)})

def sel_longitude2(data,lon1,lon2,coord='grid_xt_ocean'):
    '''Select a longitude range of data. If this longitude range loops around the -280 80 gap, roll data before selecting
    Difference from sel_longitude is that the original coordinates are restored'''
    
    if np.isnan(lon1):
        lon1 = -280
    if np.isnan(lon2):
        lon2 = 80
    if lon1<lon2:
        return data.sel({coord:slice(lon1,lon2)})
    
    n_roll = int(np.sum(data.coords[coord] >lon1))
    newdata    =  data.roll({coord:n_roll},roll_coords=True)
    
    new_coords =  np.array(newdata.coords[coord])
    new_coords[np.where(new_coords>float(lon1))]-=360
    newdata = newdata.assign_coords({coord:new_coords})
    
    out = newdata.sel({coord:slice(lon1-360,lon2)})
    old_coords = np.array(out.coords[coord])
    global_max_lon = data[coord].max().data
    old_coords[np.where(old_coords<global_max_lon-360)]+=360
    
    return out.assign_coords({coord:old_coords})

def get_bathymetry():
    return xr.open_dataset('/g/data/hh5/tmp/cosima/access-om2-01/01deg_jra55v13_iaf/output002/ocean/ocean_grid.nc').ht

def showcase(Y,Y_pred, graph = True):
    nn = ~np.isnan(Y) & ~np.isnan(Y_pred)
    MSE = np.sqrt(sklearn.metrics.mean_squared_error(Y[nn], Y_pred[nn])) #Root mean squared error
    R2 = sklearn.metrics.r2_score(Y[nn], Y_pred[nn]) #'Coefficient of determination'
    print('Root mean squared error: %.2f'%MSE)
    print('Coefficient of determination: %.2f'%R2)
    print('Pearson\'s correlation coefficient: %.2f'%np.corrcoef(Y[nn],Y_pred[nn])[1,0])

    if graph:
        # Plot outputs
        Years = np.arange(1,len(Y_pred)+1)/12
        fig = plt.figure(figsize = (10, 6)); 
        plt.scatter(Y, Y_pred,  color='black'); plt.xlabel('Test data [Sv]'); 
        plt.ylabel('Prediction [Sv]'); plt.title('Prediction vs Actual data')
        axmin = np.min((plt.ylim(),plt.xlim()))
        axmax = np.max((plt.ylim(),plt.xlim()))
        plt.xlim(axmin,axmax)
        plt.ylim(axmin,axmax)

        plt.figure(figsize=(24,6))
        plt.plot(Years,Y,c='k'); plt.plot(Years,Y_pred,c='red',linestyle='dashed'); plt.legend(['Actual','Prediction']);
        plt.xlabel('Years'); plt.ylabel('Transport [$Sv$]');
        
def evaluate_on_a_map(transect_id, data_version, data_source = '', 
                            function_str = 'RMSE', function = None, tidying = None):
    if function_str == 'RMSE':
        function = lambda A,P: np.sqrt(np.mean((np.array(A)-np.array(P))**2))
        tidying = lambda x: str(round(x/1037/10**6,2))
    elif function_str == 'R2':
        function = lambda A,B: sklearn.metrics.r2_score(A, P)
        tidying = lambda x: str(round(x,2))
    elif function_str == 'corrcoef':
        function = lambda A,P: np.corrcoef(A,P)[1,0]
        tidying = lambda x: str(round(x,2))
    elif function_str == 'id':
        def function(A,P):
            return tid
        tidying = lambda x:str(int(x))
    
    transect = xr.load_dataset('/g/data/x77/jj8842/'+data_version+'/transect_definitions.nc')
    
    metric = np.zeros(transect_id.shape)
    for i,tid in enumerate(transect_id):
        timeseries =  xr.load_dataset('/g/data/x77/jj8842/'+data_version+'/'+"{:02d}".format(int(tid))+'/'+data_source+'timeseries.nc')
        A = np.array(timeseries.iaf_actual - timeseries.iaf_actual.mean())
        P = np.array(timeseries.iaf_predicted -timeseries.iaf_predicted.mean())
        metric[i] = function(A,P)
        
    plt.figure(figsize = (20,5))
    bathymetry = get_bathymetry()
    bathymetry.sel(yt_ocean = slice(None,-25)).plot(cmap='Blues',vmin=-2000,vmax=10000)
    plot_transects(transect.sel(id=transect_id))
    plot_over_transects(metric,transect.sel(id=transect_id),tidying)
    plt.title(function_str+' '+data_source)
    
    
def fit_regression(A,B, W = None, R = None):
    '''Fit a linear regression to the dataset. 
    ie solve Ax=B and return x,
    in case the regularisation reduces the amplitude of signal
    
    A = obs/training data with dimensions of (time, gridpoints) 
    B = what you want to predict. dimension (time)
    W = weights for each timepoint. Defaults to the identity CURRENTLY NOT USED, would be ATA = A.transpose() @ W @ A and ATB = A.transpose() @ W @ B
    R = regularisation/diagonal weighting before inversion. Defaults to zeros.
        Can be a scalar or matrix (usually containing just leading diagonal elements)'''
    
    if W is None:
        W = np.ones(A.shape[0])
    if R is None:
        R = np.zeros(A.shape[1])
    if np.isscalar(R):
        R = np.ones(A.shape[1])*R
    
    
    ATA = A.transpose() @ A
    #Add regularisation
    ATA[np.diag_indices_from(ATA)] += R
    
    ATB = A.transpose() @ B
    
    C = np.linalg.solve(ATA,ATB)
    
    return C

def remove_climatology(x):
    '''Exactly what it says on the box. Remove the mean seasonal cycle
    Code adapted from CAFE-60_from_courtney.ipynb'''
    clim = x.groupby('time.month').mean('time')
    anom = (x.groupby('time.month') - clim).drop('month')
    return anom

def get_pbot_and_overturning(transect_id, expt, mascon_name,
                             remove_climatology_bool = True,
                             data_version = 'data_v1',
                             save_loc = '/g/data/x77/jj8842/',
                             freq = '1 monthly',
                             NS_extent = (-5,5),
                             ocean_cutoff = None,
                             overturning = '_overturning.nc', # Could also be some of the depth things I've generated
                             overturning_var = 'ty_trans_rho'
                            ):
    '''Returns the X and Y data required for regression or for testing effectiveness
    ocean_cutoff can be between 0 and 1; will remove mascons with an ocean proportion below ocean_cutoff'''
    
    
    #Read in the approximate mascon locations
    mascon_centres = xr.load_dataset(save_loc+'mascon_definitions/centres_'+mascon_name+'.nc')
    lons = mascon_centres.lons
    lats = mascon_centres.lats
    mascons = mascon_centres.mc
    #Shift everything onto the grid transects are defined on
    global_max_lon = 80
    lons[lons>global_max_lon] = lons[lons>global_max_lon]-360
    global_min_lon = -280
    lons[lons<global_min_lon] = lons[lons<global_min_lon]+360
    
    #Find what data we will use to predict overturning
    transects = xr.load_dataset(save_loc+data_version+'/transect_definitions.nc')
    lon_min, lon_max = transects.sel(id=transect_id).lon[0],transects.sel(id=transect_id).lon[1]
    lat_min, lat_max = transects.sel(id=transect_id).lat[0]+NS_extent[0], transects.sel(id=transect_id).lat[0]+NS_extent[1]
    
    #Find overturning
    folder = save_loc+data_version+"/{:02d}".format(int(transect_id))+'/'
    
    overturning = xr.load_dataset(folder+expt+overturning)[overturning_var]
    
    #Find pbot/mascons
    mascon = xr.load_dataset(save_loc+'regridded_pbot/'+expt+'__'+mascon_name+'.nc').sortby('mc')
    if np.isnan(lon_min):
        lon_min = global_min_lon
    if np.isnan(lon_max):
        lon_max = global_max_lon
    if lon_max>lon_min:
        bi = (lats>lat_min)&(lats<=lat_max) & (lons>lon_min) & (lons<=lon_max)
    else:
        bi = (lats>lat_min)&(lats<=lat_max) & ((lons>lon_min) | (lons<=lon_max))
        
    mascon = mascon.where(bi,drop=True)
    
    #Remove mascons which contain land
    if not (ocean_cutoff is None):
        mc_to_check = mascon_centres.sel(mc=mascon.mc)
        sufficiently_ocean = mc_to_check.where(mc_to_check.percent_ocean>ocean_cutoff,drop=True).mc
        mascon=mascon.sel(mc=sufficiently_ocean)
    elif not 'ANU' in mascon_name:
        print('You\'ve left the land in. Was this deliberate? Increase ocean_cutoff to remove mascons with less ocean proportion than that')
    
    #Check these two are comparable before I strip the time dimension
    time = overturning.time
    mascon = mascon.sel(time = time)
    assert np.all(time.values == mascon.time.values)
    
    #Strip climatology, turn everything into np arrays, haul into memory
    if remove_climatology_bool:
        Y = np.array(remove_climatology(overturning).astype('float32'))
        X = np.array(remove_climatology(mascon.ewh).astype('float32').transpose('time','mc'))
    else:
        Y = np.array((overturning-overturning.mean()).astype('float32'))
        Y = Y-np.mean(Y) #For precision reasons (I suspect it doesn't matter too much)
        X = np.array((mascon.ewh-mascon.ewh.mean('time')).astype('float32').transpose('time','mc'))
        X = X-np.mean(X,0) #Signal is so much smaller than the mean that you need to pull it off twice to maintain precision

    # Turns out the following line can definitely hurt ;( and lost me a week of confusion
    #X[:,np.isnan(X[0])] = 0 #Replace anything over land with zeros (probably no longer necessary, but can't hurt)
    not_nan = np.where(~np.isnan(X[0]))[0]
    mascon = mascon.isel(mc=not_nan)
    X = X[:,~np.isnan(X[0])]
   
    return X,Y, mascon

def noisy_R2(X,Y,C,sigma):
    Y_pred = C@X.T
    RSS = np.sum((Y_pred-Y)**2)
    NSS = np.sum((C)**2)*sigma**2*Y_pred.shape[0]
    TSS = np.sum(Y**2)
    
    return 1- (RSS+NSS)/TSS

def noisy_RMSE(X,Y,C,sigma):
    Y_pred = C@X.T
    RSS = np.sum((Y_pred-Y)**2)
    NSS = np.sum((C)**2)*sigma**2*Y_pred.shape[0]
    
    return np.sqrt((RSS+NSS)/Y_pred.shape[0])

def noisy_filtered_R2(X,Y,C,sigma,filt,filt_array):
    assert np.abs(np.sum(filt_array)-1)<0.00001
    Y_pred = C@X.T
    RSS = np.sum((filt(Y_pred)-filt(Y))**2)
    NSS = np.sum((C)**2)*sigma**2*(Y_pred.shape[0]-filt_array.shape[0])*np.sum(filt_array**2)
    TSS = np.sum(filt(Y)**2)
    
    return 1- (RSS+NSS)/TSS

def noisy_filtered_RMSE(X,Y,C,sigma,filt,filt_array):
    assert np.abs(np.sum(filt_array)-1)<0.00001
    Y_pred = C@X.T
    RSS = np.sum((filt(Y_pred)-filt(Y))**2)
    NSS = np.sum((C)**2)*sigma**2*(Y_pred.shape[0]-filt_array.shape[0])*np.sum(filt_array**2)
    TSS = np.sum(filt(Y)**2)
    
    return np.sqrt(1/Y_pred.shape[0]* (RSS+NSS))

def automated_linear_regression(namestring,transect_ids, mascon_name,
                                expt = '01deg_jra55v13_ryf9091',
                                remove_climatology_bool = True,
                                ocean_cutoff=None,
                                regularisation = 'autofit',
                                expt_test = '01deg_jra55v140_iaf_cycle3',
                                remove_climatology_test = False,
                                save_loc = '/g/data/x77/jj8842/',
                                data_version = 'data_v1',
                                freq = '1 monthly',
                                overturning = '_overturning.nc', # Could also be some of the depth things I've generated
                                overturning_var = 'ty_trans_rho',
                                NS_extent = (-5,5),
                                modify_x = None,
                                modify_y = None,
                                test_modify_x = None,
                                test_modify_y=None,
                                test_score = lambda x,y,c: sklearn.metrics.r2_score(y,c@x.T),
                                revive = False, #Whether or not to try and undo the damage done by regularisation with a constant factor
                                description = ''):
                    
    '''Given a bunch of keywords and things, find some weights that best represent how we should sum 
    mascon/bottom pressure data to get AABW transport'''
    
    for transect_id in transect_ids:
        folder = save_loc+data_version+"/{:02d}".format(int(transect_id))+'/'

        #Import data we want
        X,Y, mc = get_pbot_and_overturning(transect_id, expt, mascon_name, remove_climatology_bool, 
                                       data_version, save_loc, freq, NS_extent =NS_extent,
                                           ocean_cutoff=ocean_cutoff, overturning=overturning, overturning_var = overturning_var)


        #Adjust for sea level rise or whatever
        if not (modify_x is None):
            newX = []
            newY = []
            for f in modify_x:
                newX.append(f(X))
            for f in modify_y:
                newY.append(f(Y))
            X = np.concatenate(newX)
            Y = np.concatenate(newY)
                #If there's no data the weights are zero (avoids divide by zero erros later)

        if regularisation == 'autofit':

            X_test, Y_test, *_ = get_pbot_and_overturning(transect_id, expt_test, mascon_name, remove_climatology_bool, 
                                       data_version, save_loc, freq, NS_extent = NS_extent,
                                        ocean_cutoff=ocean_cutoff,overturning=overturning,overturning_var = overturning_var)
            if not (test_modify_x is None):
                newX = []
                newY = []
                for f in test_modify_x:
                    newX.append(f(X_test))
                for f in test_modify_y:
                    newY.append(f(Y_test))
                X_test = np.concatenate(newX)
                Y_test = np.concatenate(newY)
               
            #If either dataset is empty, then just write a bunch of zeros out as weights
            if (np.sum(np.abs(Y)) == 0) or (np.sum(np.abs(Y_test)) == 0):
                output = xr.Dataset({'weights':xr.DataArray(X[0,:]*0,dims=('mc'),coords={'mc':mc.mc})},
                                   attrs= {'expt':expt,
                                           'freq':freq,
                                           'mascon_name':mascon_name,
                                           'regularisation':regularisation,
                                           'test_expt':expt_test,
                                           'NS_extent':NS_extent,
                                           'description':description,
                                          })
                print('Writing an empty weights file for transect id: ',transect_id)
                
                output.to_netcdf(folder+namestring+'weights.nc')
                continue

            #Set up to log iterations of different regularisation
            filename = folder+namestring+'linear_fits_log.txt'
            outfile = open(filename, 'w')
            outfile.write('Regularisation, training R2, testing R2\n')

            # For different regularisation values, fit regression & evaluate
            best = [-1000,-1000,[]]
            initial_a = np.arange(-3,3,1)

            for i,a in enumerate(initial_a):
                C = fit_regression(X,Y,R=10.0**a)
                if revive:
                    m = np.std(Y)/np.std(C@X.T)
                    C = C*m

                trainR2 = test_score(X,Y,C)
                testR2  = test_score(X_test,Y_test,C)
                outfile.write(str(round(10.0**a,4))+', '+str(round(trainR2,3))+', '+str(round(testR2,3))+'\n')


                if testR2>best[0]: #If better, update
                    best[0] = testR2
                    best[1] = a
                    best[2] = C

            #For different regularisation values, fit regression & evaluate
            secondary_a = np.arange(best[1]-1,best[1]+1,0.1)
            for i,a in enumerate(secondary_a):
                C = fit_regression(X,Y,R=10.0**a)
                if revive:
                    m = np.std(Y)/np.std(C@X.T)
                    C = C*m

                trainR2 = test_score(X,Y,C)
                testR2  = test_score(X_test,Y_test,C)

                outfile.write(str(round(10.0**a,4))+', '+str(round(trainR2,3))+', '+str(round(testR2,3))+'\n')

                if testR2>best[0]: #If better, update
                    best[0] = testR2
                    best[1] = a
                    best[2] = C

            outfile.close()

            C = best[2] #Save best result for future use

                #plot timeseries
            # showcase(Y_test,(C@X_test[:,].transpose()))
            # plt.title('Regularisation = '+str(10.0**best[1])+'R$^2$ = '+str(best[0]))
            #plt.savefig(folder+'timeseries.png')

            # xr.Dataset({'iaf_actual':Y_iaf,'iaf_predicted':(C@X_iaf[:,].transpose())}).to_netcdf(
            #     folder+'timeseries.nc')

        else: #Normal regression thing
            C = fit_regression(X,Y,R=regularisation)
            
        


        #Save output
        output = xr.Dataset({'weights':xr.DataArray(C,dims=('mc'),coords={'mc':mc.mc})},
                           attrs= {'expt':expt,
                                   'freq':freq,
                                   'mascon_name':mascon_name,
                                   'regularisation':regularisation,
                                   'test_expt':expt_test,
                                   'NS_extent':NS_extent,
                                   'description':description,
                                  })

        output.to_netcdf(folder+namestring+'weights.nc')
        
def plot_mascon_weights(transect_id, data_version, name_string = '', mascon_name = '', bathymetry = True, save_loc = '/g/data/x77/jj8842/',vlim = None):
    '''Plot the weights calculated in the linear regression so I can see the spacial distribution'''
    
    weights = xr.load_dataset(save_loc+data_version+'/'+"{:02d}".format(int(transect_id))+'/'+name_string+'weights.nc').weights/1037/10**6
    plot_mascon = np.array(xr.load_dataset(save_loc+'mascon_definitions/plot_grid__'+mascon_name+'.nc').primary_mascon,int)
    transect = xr.load_dataset(save_loc+data_version+'/transect_definitions.nc').sel(id = transect_id)
    
    full_weights = np.zeros(np.max(plot_mascon)+1)
    full_weights[:] = np.nan
    full_weights[weights.mc] = weights
    gridded_weights = full_weights[plot_mascon].T
    
    mascon_centres = xr.load_dataset(save_loc+'mascon_definitions/centres_'+mascon_name+'.nc').sel(mc=weights.mc)
    lons = np.array((np.min(mascon_centres.lons)-2,np.max(mascon_centres.lons)+2))
    lons[lons>80] = lons[lons>80]-360
    lats = np.array((np.min(mascon_centres.lats-2),np.max(mascon_centres.lats+2)))
    aspect = np.abs(np.sin(np.deg2rad((float(transect.lat[0])))))
    
    plt.figure(figsize=((lons[1]-lons[0])*aspect/2+0.5,(lats[1]-lats[0])/2+0.5))
    
    if vlim is None:
        vlim = np.max((np.max(weights),-np.min(weights)))
    c = plt.imshow(gridded_weights,origin='lower',extent = (-280,80,-80,80),cmap='RdBu_r',vmin=-vlim,vmax = vlim,
                  interpolation = 'nearest', aspect = aspect,
                  )
    plt.colorbar(c)
    plt.xlim(lons)
    plt.ylim(lats)
    
    bathymetry = get_bathymetry()
    local_bath = bathymetry.sel(xt_ocean = slice(lons[0],lons[1]),yt_ocean = slice(lats[0],lats[1]))
    plt.contour(local_bath.xt_ocean,local_bath.yt_ocean,local_bath,np.arange(0,5000,1000),colors='k',alpha=0.2)
    plt.xlabel('longitude',fontsize=25)
    plt.ylabel('latitude',fontsize=25)
    plt.tight_layout()
    return

def rollave(x,w):
    '''rolling average
    x is data you're applying it to with dimensions (time, ...)
    w is window lenght'''
    y = np.zeros(x[:-w].shape)
    for i in range(y.shape[0]):
        y[i] = np.mean(x[i:i+w],0)
    return y

def low_pass(x,w):
    ''' roll a gaussian over everything
    x is data you're applying it to
    w is window length in months'''
    x_type = type(x)
    weight = xr.DataArray(np.exp(-np.linspace(-3,3,w)**2/2)/(np.sqrt(2)*np.pi), dims=["window"])
    if x_type == np.ndarray:
        x = xr.DataArray(x)
        x = x.rename({'dim_0':'time'})
    rescale = float(weight.sum('window'))
    out = x.rolling(time=w,center=True).construct(time='window').dot(weight)/rescale
    out = out.where(~np.isnan(out),drop=True)
    if x_type == np.ndarray:
        out = np.array(out)
    return out

def detrend(y,x=None):
    if x is None:
        x = np.arange(len(y))
    nn = ~np.isnan(y)
    m,c,*_ = scipy.stats.linregress(x[nn],y[nn])
    return y-m*x-c