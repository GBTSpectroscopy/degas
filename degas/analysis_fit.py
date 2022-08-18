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

    from scipy.stats import spearmanr

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

            galaxy_arr = np.append(galaxy_arr, 'all')
            corr_bin_arr =  np.append(corr_bin_arr, mybin)
            quant_arr = np.append(quant_arr,myquant)
            nval_arr = np.append(nval_arr, np.sum(idx_all))
            corr_val_arr = np.append(corr_val_arr, myresult.correlation)
            pval_arr = np.append(pval_arr, myresult.pvalue)
        
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
                
                galaxy_arr = np.append(galaxy_arr, galaxy)
                corr_bin_arr =  np.append(corr_bin_arr, mybin)
                quant_arr = np.append(quant_arr,myquant)
                nval_arr = np.append(nval_arr, np.sum(idx_gal))
                corr_val_arr = np.append(corr_val_arr, myresult.correlation)
                pval_arr = np.append(pval_arr, myresult.pvalue)
        
    

    # save results    
    result_tab = Table([galaxy_arr, corr_bin_arr, quant_arr, nval_arr, corr_val_arr, pval_arr],
                       names=('galaxy','corr_bin','quant','nval','corr_val','p_val'),
                       meta={'release':release})

    if outfile:
        result_tab.write(outfile,overwrite=True)

    return(result_tab)
                
                
