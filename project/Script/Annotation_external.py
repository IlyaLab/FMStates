import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
import scipy 
import scipy.stats as ss
import statsmodels
from statsmodels import stats
from statsmodels.stats import multitest

ROOT_DIR = os.path.abspath("../")

#### Annotation states with drug consentrations
def generate_drug_dose_annotation(samples, cellType, data_matrix_annotation_file, index_col):
    #!gsutil -m cp gs://temp_gqin/Sample1_data_MCF7_drugs_CTRP2_annotation.csv ./
    #data_matrix_annotation = pd.read_csv("Sample1_data_MCF7_drugs_CTRP2_annotation.csv") 
    data_matrix_annotation = pd.read_csv(data_matrix_annotation_file, index_col = index_col)
    #data_matrix_annotation.index = data_matrix_annotation['sig_id']
    
    #!gsutil -m cp gs://temp_gqin/gr_processing_Drug_ccl_ec50.csv ./
    EC50 = pd.read_csv(ROOT_DIR + "/Dataset/gr_processing_Drug_ccl_ec50.csv", index_col = 'Unnamed: 0')
    EC50 = EC50[EC50.columns[1:len(EC50.columns)]]
    
    phe_list = []
    data_sele = data_matrix_annotation[["pert_iname",'cell_id','pert_dose','pert_dose_unit']]
    data_sele = data_sele.loc[samples]
    phe_matrix = pd.DataFrame()
    for id in data_sele.index.values:
        if (data_sele.loc[id]['pert_iname']) in EC50.index :
            if ((np.isnan(EC50.loc[data_sele.loc[id]['pert_iname']][cellType] )) != True):
                if (data_sele.loc[id]['pert_dose_unit'] == 'µM'):
                    if (data_sele.loc[id]['pert_dose'] > EC50.loc[data_sele.loc[id]['pert_iname']][cellType]):
                        phe_list.append('up')
                    else:
                        phe_list.append('down')
          
                elif (data_sele.loc[id]['pert_dose_unit'] == 'nM'):
                    if (data_sele.loc[id]['pert_dose'] > 1000 * EC50.loc[data_sele.loc[id]['pert_iname']][cellType]):
                        phe_list.append('up')
                    else:
                        phe_list.append('down')
                else:
                    phe_list.append('unknown')
            else:
                phe_list.append('unknown')
        else:
                phe_list.append('unknown')
    phe_matrix['dosage'] = phe_list
    phe_matrix.index = data_sele.index.values
    return(phe_matrix)
# In[8]:

def plot_drug_dosage_annotation(annotation_col_1, colors, cellLIne, data_matrix_annotation_file, index_col):
    states = sorted(list(set(annotation_col_1['States'])))
    dic_state_samples = {}
    dic_state_ratio = {}
    ylist = list()
    objects = states
    for state in states:
        dic_state_samples[state] = generate_drug_dose_annotation(annotation_col_1[annotation_col_1['States'] == state].index.values,cellLIne, data_matrix_annotation_file, index_col)
        dic_state_ratio[state] = dic_state_samples[state][dic_state_samples[state]['dosage'] == 'up'].shape[0]/dic_state_samples[state][dic_state_samples[state]['dosage'] == 'down'].shape[0]
        ylist.append(dic_state_ratio[state])
    return(ylist)


def drug_response_curve(drug,a,b,c,d,fit_num_param,annotation_col_1,data_matrix_annotation_file, output_dir,colors):
    import numpy as np
    data_matrix_annotation = pd.read_csv(data_matrix_annotation_file)
    #data_matrix_annotation = pd.read_csv("Sample1_data_MCF7_drugs_CTRP2_annotation.csv")    
    data_matrix_annotation.index = data_matrix_annotation['sig_id']
    dox_data = data_matrix_annotation.loc[data_matrix_annotation['pert_iname'] ==drug ]
    states = sorted(list(set(annotation_col_1['States'])))
    dic_state_con = {}

    for state in states:
        dox_S1 = dox_data.loc[dox_data['sig_id'].isin(annotation_col_1.loc[annotation_col_1['States'] == state].index.tolist())]
        xlist_s1 = []
        for i in range(0,dox_S1.shape[0]):  
            if dox_S1.iloc[i]['pert_dose_unit'] == 'µM':
                xlist_s1.append( np.log2(dox_S1.iloc[i]['pert_dose']))
            elif dox_S1.iloc[i]['pert_dose_unit'] == 'nM':
                xlist_s1.append(np.log2(dox_S1.iloc[i]['pert_dose'] / 1000)) 
            dic_state_con[state] = xlist_s1


    import numpy as np
    plt.figure(num=None, figsize=(4, 4), dpi=300, facecolor='w', edgecolor='k')

    def drug_curve(x,fit_num_param):
        if fit_num_param == 3:
            y = ( d + c / (1 + np.exp(-1 * (x - a) / b )))
            
        if fit_num_param == 2:
            y = 1/(1 + np.exp(-1 * (x - a) / b ))
            
        return(y)
    ## drug response curve is based the orignal paper from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4718762/
    #print(y)
    t = np.arange(-10., 10., 0.1)

    #plt.plot(t, 0.1201 + 0.8798 / (1 + np.exp(-1 * (t - (-3.43)) / -2.43 )), 'k--' )
    #plt.plot(t, d + c / (1 + np.exp(-1 * (t - (a)) / b )), 'k--' )
    plt.plot(t, drug_curve(t,fit_num_param),'k--')

    i = 0
    for state in states:
        #print(state)
        if state in dic_state_con:
        #if(len(np.asarray(dic_state_con[state])) > 0):
            plt.plot(np.asarray(dic_state_con[state]), drug_curve(np.asarray(dic_state_con[state]),fit_num_param), 'o', color = colors[i]  )
        i = i + 1

    y = np.arange(0, drug_curve(a,fit_num_param), 0.1)
    x = np.asarray([a] * len(y))
    plt.plot(x, y, 'k--')

    x = np.arange(-10, a, 0.1)
    #y = np.asarray([d + c / (1 + np.exp(-1 * (0) / b ))] * len(x))
    y = np.asarray([drug_curve(a,fit_num_param)] * len(x))
    plt.plot(x, y, 'k--')
    plt.xlabel("log2_Concentration(µM)", fontdict=None, labelpad=None)
    plt.ylabel("Percentage of Valibility", fontdict=None, labelpad=None)
    plt.title(drug)
    plt.savefig(output_dir+'/'+drug+'_drugResponse_curve.png', dpi=300)
    plt.show()
    


def exclude_list(list1, list2):
    result_list = []
    for item in list1:
        if item not in list2:
            result_list.append(item)
    return(result_list)