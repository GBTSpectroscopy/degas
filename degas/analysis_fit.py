import numpy as np
from astropy.table import Table

def calc_correlation_coeffs(stack,degas_db, release='DR1',corr_bins = ['r25','mstar','ICO'],corr_quants=['ratio_HCN_CO','ratio_ltir_mean_HCN'], outfile=None):
    '''
    calculate correlation coefficients for each galaxy and quantity.

    TODO: ADD SOME SORT OF ERROR/CONFIDENCE INTERVAL ESTIMATE (VIA MONTE CARLO?)

    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    8/4/2022    A.A. Kepley     Original Code
    '''

    from scipy.stats import spearmanr, kendalltau

    release_idx = degas_db[release] == 1
    
    if not np.any(release_idx):
        print("No galaxies found in release :" + release)
        return

    galaxy_arr = []
    corr_bin_arr = []
    quant_arr = []
    nval_arr = []
    corr_val_arr = []
    pval_arr = []
    corr_val_kendall_arr = []
    pval_kendall_arr = []

    for mybin in corr_bins: 
        for myquant in corr_quants:
        
            # calculate correlation for all galaxies
            if (myquant+"_uplim" in stack.columns) & (myquant+"_lolim" in stack.columns):
                idx_all =  (stack['bin_type'] == mybin) & (stack[myquant+'_uplim'] == False) & (stack[myquant+'_lolim'] == False)                
            elif (myquant+"_uplim" in stack.columns):
                idx_all = (stack['bin_type'] == mybin) & (stack[myquant+'_uplim'] == False)
            elif (myquant+"_lolim" in stack.columns):
                                                    
                idx_all =  (stack['bin_type'] == mybin) & (stack[myquant+'_lolim'] == False)
            else:
                idx_all =  (stack['bin_type'] == mybin) 

            myresult = spearmanr(stack[idx_all]['bin_mean'], 
                                         stack[idx_all][myquant], 
                                         nan_policy='omit')

            myresult_kendall = kendalltau(stack[idx_all]['bin_mean'], 
                                          stack[idx_all][myquant], 
                                          nan_policy='omit')
            
            galaxy_arr = np.append(galaxy_arr, 'all')
            corr_bin_arr =  np.append(corr_bin_arr, mybin)
            quant_arr = np.append(quant_arr,myquant)
            nval_arr = np.append(nval_arr, np.sum(idx_all))
            corr_val_arr = np.append(corr_val_arr, myresult.correlation)
            pval_arr = np.append(pval_arr, myresult.pvalue)

            corr_val_kendall_arr = np.append(corr_val_kendall_arr, myresult_kendall.correlation)
            pval_kendall_arr = np.append(pval_kendall_arr, myresult_kendall.pvalue)
        
            # calculate correlation for each individual galaxy
            for galaxy in degas_db[release_idx]['NAME']:

                if (myquant+"_uplim" in stack.columns) & (myquant+"_lolim" in stack.columns):
                    idx_gal =  (stack['galaxy'] == galaxy) & (stack['bin_type'] == mybin) & (stack[myquant+'_uplim'] == False) & (stack[myquant+'_lolim'] == False)
                elif (myquant+"_uplim" in stack.columns):
                    idx_gal =  (stack['galaxy'] == galaxy) & (stack['bin_type'] == mybin) & (stack[myquant+'_uplim'] == False)
                elif (myquant+"_lolim" in stack.columns):
                    idx_gal =  (stack['galaxy'] == galaxy) & (stack['bin_type'] == mybin) & (stack[myquant+'_lolim'] == False)
                else:
                    idx_gal =  (stack['galaxy'] == galaxy) & (stack['bin_type'] == mybin) 

                myresult = spearmanr(stack[idx_gal]['bin_mean'], 
                                     stack[idx_gal][myquant], 
                                     nan_policy='omit')
                
                myresult_kendall = kendalltau(stack[idx_gal]['bin_mean'], 
                                              stack[idx_gal][myquant], 
                                              nan_policy='omit')

                galaxy_arr = np.append(galaxy_arr, galaxy)
                corr_bin_arr =  np.append(corr_bin_arr, mybin)
                quant_arr = np.append(quant_arr,myquant)
                nval_arr = np.append(nval_arr, np.sum(idx_gal))
                corr_val_arr = np.append(corr_val_arr, myresult.correlation)
                pval_arr = np.append(pval_arr, myresult.pvalue)
                
                corr_val_kendall_arr = np.append(corr_val_kendall_arr, myresult_kendall.correlation)
                pval_kendall_arr = np.append(pval_kendall_arr, myresult_kendall.pvalue)
    

    # save results    
    result_tab = Table([galaxy_arr, corr_bin_arr, quant_arr, nval_arr, corr_val_arr, pval_arr, corr_val_kendall_arr, pval_kendall_arr],
                       names=('galaxy','corr_bin','quant','nval','corr_val_spearman','p_val_spearman','corr_val_kendall','p_val_kendall'),
                       meta={'release':release})

    if outfile:
        result_tab.write(outfile,overwrite=True)

    return(result_tab)


                                
def fit_trend(stack, 
              bin_type='r25', 
              stack_col = 'ratio_HCN_CO', 
              exclude_gal = ''):
    '''

    Purpose: fit a linear regression to a subset of data from a stack

    TODO: 
    -- deal with log errors.
    -- how to save the models.
        -- I could just save the m and x to a file? and perhaps the fit result?
        -- from the lmfit documentation it looks like saving the model is a bit of a hassle and not guaranteed to load in different versions of python.
        -- do i want to have a driver function that does a bunch of these and saves to an astropy table??
    -- how to exclude galaxies from just one fit?? maybe move this up to driver code?? I may want to fit each galaxy at some point. Although that seems weird.
    


    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    8/19/2022   A.A. Kepley     Original Code

    '''

    # If I need to get more sophisticated try lmfit instead.
    #from scipy.stats import curve_fit
    from lmfit.models import LinearModel
    import matplotlib.pyplot as plt

    if (stack_col+'_uplim' in stack.columns) and (stack_col+'_lolim' in stack.columns):
        idx = (stack['bin_type'] == bin_type) & (stack[stack_col+'_uplim'] == False) & (stack[stack_col+'_lolim'] == False) & (stack['galaxy'] != exclude_gal)
    elif (stack_col+'_uplim' in stack.columns):
        idx = (stack['bin_type'] == bin_type) & (stack[stack_col+'_uplim'] == False) & (stack['galaxy'] != exclude_gal)
    elif (stack_col+'_lolim' in stack.columns):
        idx = (stack['bin_type'] == bin_type) & (stack[stack_col+'_lolim'] == False) & (stack['galaxy'] != exclude_gal)
        

    if bin_type == 'r25':
        xvals = np.array(stack['bin_mean'][idx])
    else:
        xvals = np.log10(np.array(stack['bin_mean'][idx]))

    yvals = np.log10(np.array(stack[stack_col][idx]))
    
    ## TODO: can't just take log values of errors.
    yvals_err = np.log10(np.array(stack[stack_col+"_err"][idx]))
    yvals_weights = 1/(yvals_err**2)

    mod = LinearModel(nan_policy='omit')
    params_mod = mod.guess(yvals,x=xvals)

    modresult = mod.fit(yvals,params=params_mod,x=xvals,verbose=True)

    print(modresult.fit_report())

    plt.clf()
    plt.scatter(xvals, yvals)
    plt.plot(xvals,modresult.best_fit,label='best fit')
    plt.xlabel(bin_type)
    plt.ylabel(stack_col)
    plt.legend()
    
    ## TODO: to get value of fit parameter and error
    ## modresult.params['slope'].value
    ## modresult.params['slope'].stderr

    return modresult
