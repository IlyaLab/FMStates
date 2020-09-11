import math
import pandas as pd
import numpy as np
import scipy 
import scipy.stats as ss
import statsmodels
from statsmodels import stats
from statsmodels.stats import multitest
import matplotlib.pyplot as plt

def plot_pval(pd_pVal, pd_cor, threshold,output_dir):
    pd_pVal = pd_pVal.astype(float)
    pd_cor = pd_cor.astype(float)
    num_sig = []
    effect_sig = []
    for factor in pd_pVal.columns.values:
        num = pd_cor.loc[ (pd_pVal[factor])<0.05].shape[0]
        num_sig.append(num)

        effect = pd_cor.loc[(pd_pVal[factor])<0.05][factor].median()
        effect_sig.append(effect)

    pd_num_effect = pd.DataFrame({'Num':num_sig,'Effect': effect_sig } )
    pd_num_effect.index = pd_pVal.columns.values

    new = pd_num_effect.sort_values(by = ['Num'], ascending=False)
    new = new.loc[new['Num'] > threshold]
    plt.figure(figsize=(6,6))
    
    #axes = new.plot.bar(rot=90, subplots=True)
    new['Num'].plot(kind='bar',fontsize = 12, color=['grey'])
    #axes = new.plot(kind='bar',subplots=True)

    #axes[1].legend(loc=1)  # doctest: +SKIP
    plt.tight_layout()
    plt.savefig(output_dir +"/Num.pdf", dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format="pdf",
        transparent=True, bbox_inches=None, pad_inches=0.2,
        frameon=None, metadata=None)
    
    plt.figure(figsize=(6,6))
    
    #axes = new.plot.bar(rot=90, subplots=True)
    new['Effect'].plot(kind='bar',fontsize = 12, color=['grey'])

    plt.tight_layout()
    plt.savefig(output_dir +"/Effect.pdf", dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format="pdf",
        transparent=True, bbox_inches=None, pad_inches=0.2,
        frameon=None, metadata=None)
    
    #plt.show()
    return(new)

def correlation_analysis(matrix_factor, Drug_IC50, min_IC50 = -1):
    '''
    Discription: 
    Get the correlation matrix and p-value matrix.
    
    Parameters:
    matrix_factor: dataframe for the FM-factors
    Drug_IC50: dataframe for the drug IC50 values
    min_IC50: threshold of IC50 for considering the effective drugs. In the example case, IC50 is log transformed. 
    '''
    
    import matplotlib.pyplot as plt
    import numpy as np
    from sklearn.linear_model import LinearRegression
    import numpy as np, statsmodels.api as sm
    import  scipy
    from scipy.stats import spearmanr, pearsonr
    
    pd_cor = pd.DataFrame()
    pd_pVal = pd.DataFrame()

    model = LinearRegression()
    for factor in matrix_factor.columns.values.tolist():
        temp_cor = []
        temp_pVal = []
        for drug in Drug_IC50.index.values:
            if np.nanmin(Drug_IC50.loc[drug]) < min_IC50:
                y = Drug_IC50.loc[drug].values
                x1 = np.array(matrix_factor[factor].values)
                x1 = x1[np.isnan(y)==False]
                y = y[np.isnan(y)==False]
                if(len(y) > 6):
                    rho, pval = spearmanr(x1,y)
                    temp_cor.append(rho)
                    temp_pVal.append(pval)
                else:
                    temp_cor.append('nan')
                    temp_pVal.append('nan')
            else:
                temp_cor.append('nan')
                temp_pVal.append('nan')
        pd_cor[factor] = temp_cor
        
        pd_pVal[factor] = temp_pVal
    pd_cor.index = Drug_IC50.index.values
    pd_pVal.index = Drug_IC50.index.values
    return(pd_cor,pd_pVal)

def get_sig_matrix(pd_pVal, pd_cor, FM_factor_list):
    '''
    Discription: 
    Get the correlation matrix with significant correlations.
    
    Parameters:
    pd_pVal: dataframe for the p-values
    pd_cor: dataframe for the correlations
    '''
    pd_pVal = pd_pVal.astype(float)
    pd_cor = pd_cor.astype(float)
    cor_matrix  =pd_cor 
    cor_matrix[cor_matrix.isna()] = 0
    cor_matrix = cor_matrix[FM_factor_list]


    cor_matrix_sig = pd_cor[pd_pVal < 0.05]
    cor_matrix_sig[cor_matrix_sig.isna()] = 0
    cor_matrix_sig = cor_matrix_sig[FM_factor_list]

    cor_matrix_simp = cor_matrix[(cor_matrix.T != 0).any()]
    cor_matrix_sig_simp = cor_matrix_sig[(cor_matrix_sig.T != 0).any()]
    
    return(cor_matrix_simp, cor_matrix_sig_simp, cor_matrix_sig)

def effective_drugs(DataFrame_DrugIC50, thereshold_percentage= 0.1, thereshold_sensitivity = -1, measurements = 6 ):
    '''
    Discription: This function is used for selecting drugs which show sensitivity(below thereshold_sensitivity) in at least of #thereshold_percentage in the population investigated. 
    
    Parameters:
    DataFrame_DrugIC50: IC50 value data matrix
    thereshold_percentage: Percentage of samples with IC50 for one drug lower than the thereshold_sensitivity. 
    thereshold_sensitivity: thereshold to consider a drug sensitive or not. IC50 smaller than this value is considered as sensitive.
    '''
    
    effect_drug = []
    Drug_ids = DataFrame_DrugIC50.index.tolist()
    for drug in range(0, DataFrame_DrugIC50.shape[0]):
        temp_effect = []
        temp_all = []
        temp = DataFrame_DrugIC50.iloc[drug].tolist()
        for i in range(0, len(temp)):
            if temp[i] < thereshold_sensitivity:
                temp_effect.append(temp[i])
            if math.isnan(float(temp[i])) == False:
                temp_all.append(Drug_ids[drug])
        if (len(temp_effect) >= len(temp_all) * thereshold_percentage) & (len(temp_all) > measurements):
            effect_drug.append(Drug_ids[drug])
    return(effect_drug)
                
def convert_toInt(list): 
    res = []
    # Converting integer list to string list 
    for i in list:
        res.append(int(i)) 
    return(res)

def convert(list): 
    res = []
    # Converting integer list to string list 
    for i in list:
        res.append(str(i)) 
    return(res)