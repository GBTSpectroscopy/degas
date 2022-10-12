import numpy as np
from astropy.table import Table
import os
import ipdb

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


def make_fit_table(stack, 
                   bins = ['r25','mstar','ICO'],
                   columns = ['ratio_HCN_CO','ratio_ltir_mean_HCN'],
                   outfile = 'test.fits',
                   exclude_gal = [],
                   plotDir = None):
    '''
    
    Purpose: run fitting code for different parameters and save
    results to table

    stack: input stack
    bins: bins to use for x vals
    columns: columns to use for y vals
    outfile: fits table for results
    exclude_gal (optional): galaxies to exclude from fit

    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    9/1/2022   A.A. Kepley     Original Code

    '''
    
    # create results table
    t = Table(names=('bin','column',
                     'slope','slope_err','intercept','intercept_err',
                     'chisq','AIC','BIC'),
              dtype=('S20','S20',
                     float, float,float,float,
                     float,float,float),
              meta=stack.meta)

    # remove any galaxies we are excluding
    if len(exclude_gal) > 0:
        for bad_gal in exclude_gal:
            idx = stack['galaxy'] != bad_gal
            stack = stack[idx]
        t.meta['EXCLUDE_GAL'] = exclude_gal

    # do fit for each bin/column combination
    for mybin in bins:

        # select out the bins we want
        idx = stack['bin_type'] == mybin
        if len(stack[idx]) == 0:
            print("Bin " + mybin+ " not found in stack. Skipping")
            continue

        # iterate over columns    
        for mycol in columns: 
            if not mycol in stack.columns:
                print("Column "+mycol+" not in stack. Skipping")
                continue
            myresult = fit_trend(stack[idx],bin_type=mybin,stack_col=mycol, plotDir=plotDir)

            t.add_row([mybin,mycol,
                       myresult.params['slope'].value,
                       myresult.params['slope'].stderr,
                       myresult.params['intercept'].value,
                       myresult.params['intercept'].stderr,
                       myresult.chisqr,
                       myresult.aic,
                       myresult.bic])

    # write out results
    t.write(outfile,overwrite=True)

    return t
    

    
                                
def fit_trend(stack, 
              bin_type='r25', 
              stack_col = 'ratio_HCN_CO',
              plotDir = None):
    '''

    Purpose: fit a linear regression to a subset of data from a stack

    TODO: 
    -- save plots for diagnostic purposes??


    Date        Programmer      Description of Changes
    ----------------------------------------------------------------------
    8/19/2022   A.A. Kepley     Original Code

    '''

    # If I need to get more sophisticated try lmfit instead.
    #from scipy.stats import curve_fit
    from lmfit.models import LinearModel
    import matplotlib.pyplot as plt

    if (stack_col+'_uplim' in stack.columns) and (stack_col+'_lolim' in stack.columns):
        idx = (stack['bin_type'] == bin_type) & (stack[stack_col+'_uplim'] == False) & (stack[stack_col+'_lolim'] == False) 
    elif (stack_col+'_uplim' in stack.columns):
        idx = (stack['bin_type'] == bin_type) & (stack[stack_col+'_uplim'] == False)  
    elif (stack_col+'_lolim' in stack.columns):
        idx = (stack['bin_type'] == bin_type) & (stack[stack_col+'_lolim'] == False) 
        

    if bin_type == 'r25':
        xvals = np.array(stack['bin_mean'][idx])
    else:
        xvals = np.log10(np.array(stack['bin_mean'][idx]))

    yvals = np.log10(np.array(stack[stack_col][idx]))
    
    ## percent errors = errors in log space
    yvals_err = stack[stack_col+"_err"][idx]/stack[stack_col][idx]
    yvals_weights = 1/(yvals_err**2)

    mygals = stack['galaxy'][idx]

    mod = LinearModel(nan_policy='omit')
    params_mod = mod.guess(yvals,x=xvals)

    #modresult = mod.fit(yvals,params=params_mod,x=xvals,verbose=True,
    #   weights = yvals_weights)

    modresult = mod.fit(yvals,params=params_mod,x=xvals,verbose=True)

    print("\n\n")
    print("**************************************************")
    print("Fit results for " + bin_type + " vs. " + stack_col)
    print(modresult.fit_report())

    plt.clf()

    gallist = np.unique(mygals)
    for gal in gallist:
        plt.errorbar(xvals[gal == mygals], yvals[gal == mygals],
                     yerr = yvals_err[gal==mygals],
                     label=gal)

    plt.plot(xvals,modresult.best_fit,label='best fit',color='darkgray',linewidth=3.0)
    plt.xlabel(bin_type)
    plt.ylabel(stack_col)
    plt.legend()
        
    if plotDir:
        if not os.path.exists(plotDir):
            os.mkdir(plotDir)
        plt.savefig(os.path.join(plotDir,bin_type + "_" + stack_col + "_fit.png"))
        plt.savefig(os.path.join(plotDir,bin_type + "_" + stack_col + "_fit.pdf"))

    return modresult
