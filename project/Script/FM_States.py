## load libriaries

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


#### User input: User can select functional moduels and input the gene expression matrix 

ROOT_DIR = os.path.abspath("../")

#### Define functional module factors. Here we defined up_regulation_strength, down_regulation_strength, ssGSEA, and TF_strength as shown below. User can also define their own factors. 

def intersection(lst1, lst2): 
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3 

##### Need to add a user defined function pattern in the "generate_factor" function. 
def ratio_up(x):
    result = len(x[x > 1.6])/len(x)
    return(result)

def factor_up(data_matrix, GS_20,KEGG_modules):
    matrix_factor_up = pd.DataFrame()
    for i in range(0,len(GS_20)):
        matrix_module = data_matrix[intersection(data_matrix.columns, list(KEGG_modules[KEGG_modules['name'] == GS_20[i]]['member']))]
        
        module_up_result = []
        #for sample in range(0,matrix_module.shape[0]):
        for sample in matrix_module.index.values:
            x = matrix_module.loc[sample]
            module_up_result.append(ratio_up(x))
        matrix_factor_up[GS_20[i]] = module_up_result
    #matrix_factor_up.index.names  = data_matrix_MCF7_CTRP2.index
    
    matrix_factor_up.columns = matrix_factor_up.columns.values + '_up'
    matrix_factor_up.index = list(data_matrix.index.values)
    return(matrix_factor_up)

def factor_up_absolute(data_matrix, GS_20,KEGG_modules):
    matrix_factor_up = pd.DataFrame()
    data_matrix_zscore = ss.zscore(data_matrix, axis=1)
    #data_matrix_zscore = ss.zscore(data_matrix_zscore, axis=0)
    data_matrix_zscore_df = pd.DataFrame(data=data_matrix_zscore,    # values
                                        index=list(data_matrix.index.values),    # 1st column as index
                                      columns=data_matrix.columns.values )
    
    for i in range(0,len(GS_20)):
        matrix_module = data_matrix_zscore_df[intersection(data_matrix_zscore_df.columns, list(KEGG_modules[KEGG_modules['name'] == GS_20[i]]['member']))]
        #matrix_module.index = data_matrix_MCF7_CTRP2['Unnamed: 0']
        
        module_up_result = []
        #for sample in range(0,matrix_module.shape[0]):
        for sample in matrix_module.index.values:
            x = matrix_module.loc[sample]
            module_up_result.append(ratio_up(x))
        matrix_factor_up[GS_20[i]] = module_up_result
    #matrix_factor_up.index.names  = data_matrix_MCF7_CTRP2.index
    
    matrix_factor_up.columns = matrix_factor_up.columns.values + '_up'
    matrix_factor_up.index = list(data_matrix.index.values)
    return(matrix_factor_up)


#Define factor 2: down_regulation_strength

def ratio_down(x):
    result = -1 * len(x[x < -1.6])/len(x)
    return(result)

def factor_down(data_matrix, GS_20,KEGG_modules):
    matrix_factor_down = pd.DataFrame()
    for i in range(0,len(GS_20)):
        matrix_module = data_matrix[intersection(data_matrix.columns,list(KEGG_modules[KEGG_modules['name'] == GS_20[i]]['member']))]
        #matrix_module.index = data_matrix_MCF7_CTRP2['Unnamed: 0']
        
        module_down_result = []
        #for sample in range(0,matrix_module.shape[0]):
        for sample in matrix_module.index.values:
            x = matrix_module.loc[sample]
            module_down_result.append(ratio_down(x))
        
        matrix_factor_down[GS_20[i]] = module_down_result

    matrix_factor_down.columns = matrix_factor_down.columns.values+'_down'
    matrix_factor_down.index = list(data_matrix.index.values)
    #print(matrix_factor_down.index.values)
    return(matrix_factor_down)

def factor_down_absolute(data_matrix, GS_20,KEGG_modules):
    data_matrix_zscore = ss.zscore(data_matrix, axis=1)
    #data_matrix_zscore = ss.zscore(data_matrix_zscore, axis=0)
    data_matrix_zscore_df = pd.DataFrame(data=data_matrix_zscore,    # values
                                        index=list(data_matrix.index.values),    # 1st column as index
                                      columns=data_matrix.columns.values )
    
    matrix_factor_down = pd.DataFrame()
    for i in range(0,len(GS_20)):
        matrix_module = data_matrix_zscore_df[intersection(data_matrix_zscore_df.columns,list(KEGG_modules[KEGG_modules['name'] == GS_20[i]]['member']))]
        #matrix_module.index = data_matrix_MCF7_CTRP2['Unnamed: 0']
        
        module_down_result = []
        #for sample in range(0,matrix_module.shape[0]):
        for sample in matrix_module.index.values:
            x = matrix_module.loc[sample]
            module_down_result.append(ratio_down(x))
        
        matrix_factor_down[GS_20[i]] = module_down_result

    matrix_factor_down.columns = matrix_factor_down.columns.values+'_down'
    matrix_factor_down.index = list(data_matrix.index.values)
    #print(matrix_factor_down.index.values)
    return(matrix_factor_down)

#Define factor 3: single sample GSEA enrichment score

def factor_ssGSEA(data_matrix, GS_20,KEGG_modules):
    
    from GSVA import gsva, gmt_to_dataframe
    import math
    limit = 1000
    n = math.ceil(data_matrix.shape[0] / limit)
    
    module_selected_gmt = KEGG_modules.loc[KEGG_modules['name'].isin( GS_20) ]
    
    pathways_df = gsva(data_matrix.T,module_selected_gmt,mx_diff = False, parallel_sz = 2, method = 'ssgsea')
    
    matrix_factor_gsva = pd.DataFrame()
    result_list = []
    for i in range(0,n):
        print(i)
        data = data_matrix.iloc[i*1000: i*1000+1000]
        pathways_df = gsva(data.T,module_selected_gmt,mx_diff = False, parallel_sz = 2, method = 'ssgsea')
        result = pathways_df.T
        result.columns = result.columns.values + '_ssGSEA'
        result.index = list(data.index.values)
        result_list.append(result)
    matrix_factor_gsva = pd.concat(result_list, axis = 0)

    #pathways_df.T.columns = pathways_df.T.columns.values+'_ssGSEA'
    #result = pathways_df.T
    #result.columns = result.columns.values + '_ssGSEA'
    #result.index = list(data_matrix.index.values)
    #print(result.index.values)
    return(matrix_factor_gsva)

#Define factor 4: transcriptional strength score

def factor_tf(data_matrix_MCF7_CTRP2,GS_20,TFs): #remove factor activation = True
    matrix_factor_tf = pd.DataFrame()
    #TFs = pd.read_csv("/Users/gqin/Documents/ISB/FM_CStates/Data/TF_fromGR_enriched.csv")
    for module in range(0,len(GS_20)):
        expr_tf = data_matrix_MCF7_CTRP2[intersection(list(set(TFs.loc[TFs['Pathway'] == GS_20[module]]['TFs'])),data_matrix_MCF7_CTRP2.columns.values )]
        tf_module = intersection(list(set(TFs.loc[TFs['Pathway'] == GS_20[module]]['TFs'])),expr_tf.columns.values )
        TFs_g = TFs.loc[TFs['TFs'].isin(tf_module)]
        TFs_g_m = TFs_g.loc[TFs_g['Pathway'] == GS_20[module]]

        expr_tf_order = expr_tf[list(TFs_g_m['TFs'])]
        f = TFs_g_m['Score']/sum(TFs_g_m['Score'])

        tf_list = list()

        for sample in range(0,expr_tf_order.shape[0]):
            sum_tf = 0
            for gene in range(0,len(expr_tf_order.iloc[1])):
                sum_tf = sum_tf + expr_tf_order.iloc[sample][gene] * f.iloc[gene]
            tf_list.append(sum_tf)
           
        matrix_factor_tf[GS_20[module]] = tf_list
    
    #if activation  == True:
    #    matrix_factor_tf.columns = matrix_factor_tf.columns.values + '_tf_act'
    #elif activation == False:
    #    matrix_factor_tf.columns = matrix_factor_tf.columns.values + '_tf_sup'
    
    matrix_factor_tf.columns = matrix_factor_tf.columns.values + '_tf'
    
    matrix_factor_tf.index = data_matrix_MCF7_CTRP2.index
    return(matrix_factor_tf)


def generate_factor(data_matrix_MCF7_CTRP2, GS_20, dic_module, TFs, UP = False, DOWN = False, ssGSEA = False, TF = False, absolute = False): ##TFs_active, TFs_supress,
    factor_matrix_list = list()
    matrix_factor_up = pd.DataFrame()
    matrix_factor_down = pd.DataFrame()
    matrix_factor_ssGSEA = pd.DataFrame()
    matrix_factor_tf = pd.DataFrame()
    if UP == True:
        if absolute == False:
            matrix_factor_up = factor_up(data_matrix_MCF7_CTRP2, GS_20, dic_module)
        elif absolute == True:
            matrix_factor_up = factor_up_absolute(data_matrix_MCF7_CTRP2, GS_20, dic_module)
        factor_matrix_list.append(matrix_factor_up)
        
    if DOWN == True:
        if absolute == False:
            matrix_factor_down = factor_down(data_matrix_MCF7_CTRP2, GS_20, dic_module)
        elif absolute == True:
            matrix_factor_down = factor_down_absolute(data_matrix_MCF7_CTRP2, GS_20, dic_module)
        factor_matrix_list.append(matrix_factor_down)
        
    if ssGSEA == True:
        matrix_factor_ssGSEA = factor_ssGSEA(data_matrix_MCF7_CTRP2,  GS_20, dic_module)
        factor_matrix_list.append(matrix_factor_ssGSEA)
        
    if TF == True:
        
        matrix_factor_tf = factor_tf(data_matrix_MCF7_CTRP2, GS_20, TFs)
        factor_matrix_list.append(matrix_factor_tf)
        
        #matrix_factor_tf_active = factor_tf(data_matrix_MCF7_CTRP2, GS_20, TFs_active, activation = True)
        #factor_matrix_list.append(matrix_factor_tf_active)
        
        #matrix_factor_tf_supress = factor_tf(data_matrix_MCF7_CTRP2, GS_20, TFs_supress, activation = False)
        #factor_matrix_list.append(matrix_factor_tf_supress)
    if len(factor_matrix_list) == 0:
        print("Error!! No factor is selected!")
        
    #matrix_factor = pd.concat([matrix_factor_up, matrix_factor_down, matrix_factor_ssGSEA, matrix_factor_tf], axis=1)
    matrix_factor = pd.concat(factor_matrix_list, axis = 1)
    return(matrix_factor)

def Rank_bySample(x):
    x_order = ss.rankdata(x)/len(x)*100
    #result = len(x[x > (1 - 0.05)*100])
    return(x_order)

def ratio_up_percentile(x):
    result = len(x[x > (1 - 0.05)*100])/len(x)
    return(result)

def ratio_down_percentile(x):
    result = -1 * len(x[x < 0.05 *100])/len(x)
    return(result)

def factor_down_percentile(matrix_ranked, GS_20,KEGG_modules):
    matrix_factor_down = pd.DataFrame()
    for i in range(0,len(GS_20)):
        matrix_module = matrix_ranked[intersection(matrix_ranked.columns, list(KEGG_modules[KEGG_modules['name'] == GS_20[i]]['member']))]
        #matrix_module.index = data_matrix_MCF7_CTRP2['Unnamed: 0']

        module_down_result = []
        #for sample in range(0,matrix_module.shape[0]):
        for sample in matrix_module.index.values:
            x = matrix_module.loc[sample]
            module_down_result.append(ratio_down_percentile(x))
        matrix_factor_down[GS_20[i]] = module_down_result
        #matrix_factor_up.index.names  = data_matrix_MCF7_CTRP2.index
        
    matrix_factor_down.columns = matrix_factor_down.columns.values + '_down'
    matrix_factor_down.index = list(matrix_ranked.index.values)
    return(matrix_factor_down)

def factor_up_percentile(matrix_ranked, GS_20,KEGG_modules):
    matrix_factor_up = pd.DataFrame()
    for i in range(0,len(GS_20)):
        matrix_module = matrix_ranked[intersection(matrix_ranked.columns, list(KEGG_modules[KEGG_modules['name'] == GS_20[i]]['member']))]
        #matrix_module.index = data_matrix_MCF7_CTRP2['Unnamed: 0']

        module_up_result = []
        #for sample in range(0,matrix_module.shape[0]):
        for sample in matrix_module.index.values:
            x = matrix_module.loc[sample]
            module_up_result.append(ratio_up_percentile(x))
        matrix_factor_up[GS_20[i]] = module_up_result
        #matrix_factor_up.index.names  = data_matrix_MCF7_CTRP2.index
        
    matrix_factor_up.columns = matrix_factor_up.columns.values + '_up'
    matrix_factor_up.index = list(matrix_ranked.index.values)
    return(matrix_factor_up)

#### Measuring factor difference between two groups

def cohen_dist(vec1, vec2):
    M1 = np.mean(vec1)
    N1 = len(vec1)
    M2 = np.mean(vec2)
    N2 = len(vec2)
    SD_pooled = np.sqrt(
                        (np.square(np.std(vec1)) * (N1 - 1) 
                        + np.square(np.std(vec2)) * (N2 - 1)) 
                        / (N1 + N2 -2)
                       )
    d = (M1 - M2) / SD_pooled
    return(d)


#Compare whether one factor show significant difference between two groups. 
def Diff_comp(ctrl_factor,matrix_factor):
    p_list = []
    d_list = []
    module_list = []
    for module in matrix_factor.columns.values:
        if sum(ctrl_factor[module]) != 0 :
            module_list.append(module)
            static_p = scipy.stats.mannwhitneyu(ctrl_factor[module], matrix_factor[module], use_continuity=True, alternative=None)
            p = static_p[1]
            p_list.append(p)
            M1 = np.mean(matrix_factor[module])
            N1 = matrix_factor.shape[0]
            M2 = np.mean(ctrl_factor[module])
            N2 = ctrl_factor.shape[0]
            SD_pooled = np.sqrt(
                (np.square(np.std(matrix_factor[module])) * (N1-1)
                 + np.square(np.std(ctrl_factor[module])) * (N2-1)) 
                / (N1 + N2 - 2)
                )
            d = (M1 - M2)/SD_pooled
            d_list.append(d)
        
    FDR_result = statsmodels.stats.multitest.multipletests(p_list, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    FDR_list = FDR_result[1]
    FDR_list = -1 * np.log10(FDR_list)
    diff_matrix = {'SizeEffect':d_list, 'P':p_list, '-log10(FDR)':FDR_list} 
    diff_dataframe = pd.DataFrame(diff_matrix)
    diff_dataframe.index = module_list
    return(diff_dataframe)


#### Normalize functional module factors using the reference factors
def compare_factors_between_experimental_to_reference(ctrl_factor, exp_factor):
    '''
    Compare factors between the experimental group with the reference group.
    '''
    compare_factor_matrix = exp_factor
    for factor in ctrl_factor.columns:
        if '_up' in factor:
            x = ctrl_factor[factor]
            num_item = len(x)
            arr = []
            for i in range(0,len(exp_factor[factor])):
                item = (exp_factor[factor][i])
                ref = x.values
                y = np.append(ref,item)
                order = y.argsort()
                ranks = order.argsort()    
                new_p = ranks[-1]/(num_item)
                arr.append(new_p)
            compare_factor_matrix[factor] = arr
                
        if '_down' in factor:
            x =  ctrl_factor[factor]
            num_item = len(x)
            arr = []
            for i in range(0,len(exp_factor[factor])):
                item = (exp_factor[factor][i])
                ref = x.values
                y = np.append(ref,item)
                order = y.argsort()
                ranks = order.argsort()    
                new_p = (num_item - ranks[-1])/num_item
                arr.append(new_p)
            compare_factor_matrix[factor] = arr
            
        if ( '_ssGSEA'  in factor ) or ('_tf' in factor):
            compare_factor_matrix[factor] = (compare_factor_matrix[factor] - np.mean(ctrl_factor[factor]))/np.std(ctrl_factor[factor])

    return(compare_factor_matrix)

#### Annotation different states
def Get_features_one_vs_all_others(data, 
                                   threshold, 
                                   threshold_SE, 
                                   annotation_col_1, 
                                   select_modules): 
    '''
    The function of Get_features_one_vs_all_others is used to figure out which FM-factors show different across all the states.

    Parameters:
    data: "FM fator matrix"
    threshold: "Threshold for p-value, eg. 0.01"
    threshold_SE: "Threshold for the effect size, eg. 1", 
    annotation_col_1: "annotation of the states",
     elect_modules:  "modules selected for analysis"
    '''

    data['States'] = annotation_col_1.loc[data.index.values]
    sizeeffect_matrix = pd.DataFrame()
    states = list(set(annotation_col_1['States']))
    factors = []
    dic_factor = {}
    for state in states:
        factor = state + '_up'
        factors.append(factor)
        factor = state + '_down'
        factors.append(factor)
        factor = state + '_up_tf'
        factors.append(factor)
        factor = state + '_down_tf'
        factors.append(factor)

    for factor in factors:
        dic_factor[factor] = set()

    dic_comp = {}
    for Factor in select_modules:
        dic_comp[Factor] = {}
        for state in states:
            dic_comp[Factor][state] = data[data['States'] == state][Factor].values

        label = True 
        for state1 in states:
            for state2 in states:
                if state1 != state2:
                    if scipy.stats.ranksums(dic_comp[Factor][state1], dic_comp[Factor][state2])[1] > threshold:
                        label = False
                        
        if label == True:
            cohen_list_tmp = []
            for state1 in states:
                t1 = dic_comp[Factor][state1]
                t2 = []
                for state2 in states:
                    if state2 != state1:
                        for value in dic_comp[Factor][state2]:
                            t2.append(value)
                cohen_list_tmp.append((cohen_dist(t1,t2)))

            sizeeffect_matrix[Factor] = cohen_list_tmp
            sizeeffect_matrix.index = states

    labels = sizeeffect_matrix.columns
    dic_state = {}

    for state in states:
        dic_state[state] = {}
        dic_state[state]['up'] = []
        dic_state[state]['down'] = []
        dic_state[state]['down_tf'] = []
        dic_state[state]['up_tf'] = []

        for i in range(0,len(sizeeffect_matrix.loc[state])):
            if sizeeffect_matrix.loc[state][i] > threshold_SE:
                Factor = labels[i]
                items = Factor.split('_')
                if items[1] == 'ssGSEA':
                    dic_state[state]['up'].append(items[0])
                if items[1] == 'up':
                    dic_state[state]['up'].append(items[0])
                if items[1] == 'down':
                    dic_state[state]['down'].append(items[0])
                if items[1] == 'tf':
                    dic_state[state]['up_tf'].append(items[0])

            if sizeeffect_matrix.loc[state][i] < -1 * threshold_SE:
                Factor = labels[i]
                items = Factor.split('_')
                if items[1] == 'ssGSEA':
                    dic_state[state]['down'].append(items[0])
                if items[1] == 'up':
                    dic_state[state]['down'].append(items[0])
                if items[1] == 'down':
                    dic_state[state]['up'].append(items[0])
                if items[1] == 'tf':
                    dic_state[state]['down_tf'].append(items[0])
    for state in states:
        print(state)
        if len(dic_state[state]['up']) > 0:
            print("up regulated in : ")
            print(','.join(list(set(dic_state[state]['up']))))
        if len(dic_state[state]['up_tf'])> 0:
            print("up regulated in TF: ") 
            print(','.join(list(set(dic_state[state]['up_tf']))))
        if len(dic_state[state]['down']) > 0:   
            print("down regulated in : ")
            print(','.join(list(set(dic_state[state]['down']))))
        if len(dic_state[state]['down_tf']) > 0:   
            print("down regulated in TF: ")
            print(','.join(list(set(dic_state[state]['down_tf']))))
        print('\n') 
    return(sizeeffect_matrix)


def Get_features_one_vs_one(data, 
                 threshold, 
                 threshold_SE, 
                 threshold_num_of_states_show_difference, 
                 annotation_col_1, 
                 select_modules): 
    '''
    The function of Get_features_one_vs_one is used to figure out which FM-factors show different acros all the states.
    Parameters: 
        data: "FM fator matrix"
        threshold: "Threshold for p-value"
        threshold_SE: "Threshold for the effect size", 
        threshold_num_of_states_show_difference: "Threshold_num_of_states_show_difference",
        annotation_col_1: "annotation of the states",
        select_modules: "modules selected for analysis"
    '''
    data['States'] = annotation_col_1.loc[data.index.values]
    sizeeffect_matrix = pd.DataFrame()
    states = list(set(annotation_col_1['States']))
    factors = []
    dic_factor = {}
    for state in states:
        factor = state + '_up'
        factors.append(factor)
        factor = state + '_down'
        factors.append(factor)
        factor = state + '_up_tf'
        factors.append(factor)
        factor = state + '_down_tf'
        factors.append(factor)

    for factor in factors:
        dic_factor[factor] = set()
        
    dic_state = {}
    for state in states:
        dic_state[state] = {}
        dic_state[state]['up'] = []
        dic_state[state]['down'] = []
        dic_state[state]['down_tf'] = []
        dic_state[state]['up_tf'] = []
        
    dic_comp = {}
    for Factor in select_modules:
        dic_comp[Factor] = {}
        for state in states:
            dic_comp[Factor][state] = data[data['States'] == state][Factor].values

        label = []
        dic_sig = {}
        for state1 in states:
            count = 0
            for state2 in states:
                if state1 != state2:
                    if scipy.stats.ranksums(dic_comp[Factor][state1], dic_comp[Factor][state2])[1] > threshold:
                        label.append(False)
                    else:
                        label.append(True)
                        count = count + 1
            
            #print(count)
            if count >= threshold_num_of_states_show_difference:
                if Factor not in dic_sig:
                    dic_sig[Factor] = [state1]
                else: 
                    dic_sig[Factor].append(state1)
                    
                cohen_list_tmp = []
                t1 = dic_comp[Factor][state1]
                t2 = []
                for state2 in states:
                    if state2 != state1:
                        for value in dic_comp[Factor][state2]:
                            t2.append(value)
                    cohen_list_tmp.append((cohen_dist(t1,t2)))
                
                
                se_label = False
                if cohen_dist(t1,t2) > threshold_SE:
                    items = Factor.split('_')
                    se_label = True
                    if items[1] == 'ssGSEA':
                        dic_state[state1]['up'].append(items[0])
                    if items[1] == 'up':
                        dic_state[state1]['up'].append(items[0])
                    if items[1] == 'down':
                        dic_state[state1]['down'].append(items[0])
                    if items[1] == 'tf':
                        dic_state[state1]['up_tf'].append(items[0])
        
                    
                elif cohen_dist(t1,t2) < -1 * threshold_SE:
                    items = Factor.split('_')
                    se_label = True
                    if items[1] == 'ssGSEA':
                        dic_state[state1]['down'].append(items[0])
                    if items[1] == 'up':
                        dic_state[state1]['down'].append(items[0])
                    if items[1] == 'down':
                        dic_state[state1]['up'].append(items[0])
                    if items[1] == 'tf':
                        dic_state[state1]['down_tf'].append(items[0])
                 
                if se_label == True:
                    cohen_list_tmp = []
                    for state1 in states:
                        t1 = dic_comp[Factor][state1]
                        t2 = []
                        for state2 in states:
                            if state2 != state1:
                                for value in dic_comp[Factor][state2]:
                                    t2.append(value)
                        cohen_list_tmp.append((cohen_dist(t1,t2)))

                    sizeeffect_matrix[Factor] = cohen_list_tmp
                    sizeeffect_matrix.index = states
                
    for state in states:
        print(state)
        if len(dic_state[state]['up']) > 0:
            print("up regulated in : ")
            print(','.join(list(set(dic_state[state]['up']))))
        if len(dic_state[state]['up_tf'])> 0:
            print("up regulated in TF: ") 
            print(','.join(list(set(dic_state[state]['up_tf']))))
        if len(dic_state[state]['down']) > 0:   
            print("down regulated in : ")
            print(','.join(list(set(dic_state[state]['down']))))
        if len(dic_state[state]['down_tf']) > 0:   
            print("down regulated in TF: ")
            print(','.join(list(set(dic_state[state]['down_tf']))))
        print('\n') 
    return([sizeeffect_matrix,dic_state])


def get_tfpairs_for_select_pathways_directed(data_matrix_MCF7_CTRP2,pathway_select,dic_module, active = True): ## remove parameter active = True
    
    Resulting_pairs = get_tf_target_pair(data_matrix_MCF7_CTRP2,pathway_select,dic_module)
    if active == True:
        Resulting_pairs = Resulting_pairs.loc[Resulting_pairs['corr'] > 0]
    elif active == False:
        Resulting_pairs = Resulting_pairs.loc[Resulting_pairs['corr'] < 0]
        
    result_df = pd.DataFrame()
    
    for i in range(0,len(pathway_select)):
        pathway = pathway_select[i]
        test = Resulting_pairs.loc[Resulting_pairs['Target'].isin(dic_module[pathway])]
        pathway_list = []
        TF_list = []
        Occupancy_inPathway_list = []
        Occupancy_inTF_list = []
        Score_list = []

        for gene in list(set(test['TF'])):
            A = len(test.loc[test['TF'] == gene])
            B = test.shape[0] - A
            C = Resulting_pairs[Resulting_pairs['TF'] == gene].shape[0] - A
            D = Resulting_pairs.shape[0] - B - C + A
            oddsratio, pvalue = scipy.stats.fisher_exact([[A, B], [C, D]], alternative= 'greater')
            
            if pvalue < 0.05 and oddsratio > 0:
                Occupancy_inPathway = test.loc[test['TF'] == gene].shape[0]/test.shape[0]
                #Occupancy_inPathway estimate the percent of genes regulated by this transcription factor. The higher values means more enriched in this pathway.
                Occupancy_inTF = test.loc[test['TF'] == gene].shape[0]/Resulting_pairs.loc[Resulting_pairs['TF'] == gene].shape[0]
                #Occupancy_inTF estimate the percent of target genes for this transcription factor in this pathway. The higher values means more enriched in this pathway.
                Score = Occupancy_inPathway * Occupancy_inTF  
                pathway_list.append(pathway)
                TF_list.append(gene)
                Occupancy_inPathway_list.append(Occupancy_inPathway)
                Occupancy_inTF_list.append(Occupancy_inTF)
                Score_list.append(Score)

        result = pd.DataFrame({"Pathway":pathway_list, "TFs":TF_list, "Occupancy_inPathway":Occupancy_inPathway_list, "Occupancy_inTF":Occupancy_inTF_list, "Score":Score_list})
        result_df = pd.concat([result_df, result])
    return(result_df)

def TF_annotation(x,TF_pairs,annotation_col_1, output_dir, active = True):
    TF_used_df = TF_pairs[TF_pairs['Occupancy_inPathway'] > 0.03]
    TF_used_df = TF_used_df.loc[TF_used_df['Occupancy_inTF'] > 0.2]
    #TF_used_df = TF_pairs_active
    pathway_used = set()
    for item in x.columns.tolist():
        pathway_used.add(item.split('_')[0])
    pathway_used = list(pathway_used)
    TF_used_df = TF_used_df.loc[TF_used_df['Pathway'].isin(pathway_used)]
    TF_used = set(TF_used_df['TFs'])

    data_matrix_MCF7_CTRP2 = pd.read_csv(os.path.join(ROOT_DIR, "Data/Sample1_data_MCF7_drugs_CTRP2.csv"))
    data_matrix_MCF7_CTRP2 = data_matrix_MCF7_CTRP2.set_index('Unnamed: 0')
    TF_used = TF_used.intersection(set(data_matrix_MCF7_CTRP2.columns.tolist()))

    dic_sample = {}
    for state in list(set(annotation_col_1['States'])):
        dic_sample[state] = annotation_col_1.loc[annotation_col_1['States'].isin([state])].index.tolist()

    TF_list = list(TF_used)
    p_list = []
    for TF in TF_list:
        p_list.append(scipy.stats.f_oneway( * [data_matrix_MCF7_CTRP2[TF][dic_sample[mydata]].values for mydata in dic_sample])[1])

    pvals_corrected = statsmodels.stats.multitest.multipletests(p_list, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
    FDR_list = (pvals_corrected[1]).tolist()

    FDR_df = pd.DataFrame()

    FDR_df['TF'] = TF_list
    FDR_df['FDR'] = FDR_list
    
    sizeeffect_matrix = pd.DataFrame()
    states = list(set(annotation_col_1['States']))
    for Factor in TF_list:
        cohen_list_tmp = []
        for state1 in states:

            t1 = data_matrix_MCF7_CTRP2[Factor][dic_sample[state1]].values
            t2 = []
            for state2 in states:
                if state2 != state1:
                    for value in data_matrix_MCF7_CTRP2[Factor][dic_sample[state2]].values:
                        t2.append(value)

            cohen_list_tmp.append((cohen_dist(t1,t2)))
        sizeeffect_matrix[Factor] = cohen_list_tmp
        sizeeffect_matrix.index = states
    
    if active == True:
        sizeeffect_matrix.transpose().to_csv(output_dir + "/TF_effect_active.csv")
        TF_path_interaction = TF_used_df.loc[TF_used_df['TFs'].isin(list(FDR_df.loc[FDR_df['FDR'] < 0.05]['TF'].values))]
        TF_path_interaction.to_csv(output_dir+"/TF_path_active.csv")

        node_par = pd.DataFrame()
        node_par['id'] = list(set(TF_path_interaction['TFs'])) + list(set(TF_path_interaction['Pathway']))
        node_par['atr'] = ['TF']*len(set(TF_path_interaction['TFs'])) + ['Pathway'] * len(set(TF_path_interaction['Pathway']))
        node_par.to_csv(output_dir + "/TF_path_atrr_active.csv")
    elif active == False:
        sizeeffect_matrix.transpose().to_csv(output_dir + "/TF_effect_suppress.csv")
        TF_path_interaction = TF_used_df.loc[TF_used_df['TFs'].isin(list(FDR_df.loc[FDR_df['FDR'] < 0.05]['TF'].values))]
        TF_path_interaction.to_csv(output_dir + "/TF_path_suppress.csv")

        node_par = pd.DataFrame()
        node_par['id'] = list(set(TF_path_interaction['TFs'])) + list(set(TF_path_interaction['Pathway']))
        node_par['atr'] = ['TF']*len(set(TF_path_interaction['TFs'])) + ['Pathway'] * len(set(TF_path_interaction['Pathway']))
        node_par.to_csv(output_dir + "/TF_path_atrr_suppress.csv")
    return(TF_path_interaction)


#### Annotation states with drug consentrations
def generate_drug_dose_annotation(samples, cellType):
    #!gsutil -m cp gs://temp_gqin/Sample1_data_MCF7_drugs_CTRP2_annotation.csv ./
    data_matrix_annotation = pd.read_csv("Sample1_data_MCF7_drugs_CTRP2_annotation.csv")    
    data_matrix_annotation.index = data_matrix_annotation['sig_id']
    
    #!gsutil -m cp gs://temp_gqin/gr_processing_Drug_ccl_ec50.csv ./
    EC50 = pd.read_csv("gr_processing_Drug_ccl_ec50.csv")
    EC50.index = EC50['Unnamed: 0']
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

def plot_drug_dosage_annotation(annotation_col_1, colors, cellLIne):
    states = sorted(list(set(annotation_col_1['States'])))
    dic_state_samples = {}
    dic_state_ratio = {}
    ylist = list()
    objects = states
    for state in states:
        dic_state_samples[state] = generate_drug_dose_annotation(annotation_col_1[annotation_col_1['States'] == state].index.values,cellLIne)
        dic_state_ratio[state] = dic_state_samples[state][dic_state_samples[state]['dosage'] == 'up'].shape[0]/dic_state_samples[state][dic_state_samples[state]['dosage'] == 'down'].shape[0]
        ylist.append(dic_state_ratio[state])
    return(ylist)


def drug_response_curve(drug,a,b,c,d,fit_num_param,annotation_col_1):
    import numpy as np
    data_matrix_annotation = pd.read_csv("Sample1_data_MCF7_drugs_CTRP2_annotation.csv")    
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
        print(state)
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
    plt.xlabel("log2_IC50(µM)", fontdict=None, labelpad=None)
    plt.ylabel("Percentage of Valibility", fontdict=None, labelpad=None)
    plt.title(drug)
    plt.savefig(output_dir+'/'+drug+'_drugResponse_curve.png', dpi=300)
    plt.show()
    
def FDA_drugs_list(drug_list):
    #!gsutil -m cp gs://temp_gqin/FDA_drugs_2019.5.7.txt ./
    FDA_drugs = pd.read_csv("FDA_drugs_2019.5.7.txt", sep='\t', header = None)
    FDA_drugs.columns = ["drug"]
    IN_list = []
    OUT_list = []
    result = pd.DataFrame()
    for drug in drug_list:
        if (drug in FDA_drugs['drug'].values):
            IN_list.append("In")
        else:
            turn = "OFF"
            for id in FDA_drugs['drug'].values:
                #print(id)
                if drug in id: 
                    turn = "ON"
                    break
            if turn == "ON":
                IN_list.append("In") 
            else:
                IN_list.append("Out") 
    result = pd.DataFrame(IN_list)  
    result.columns = ["FDA_approved"]
    result.index = drug_list
    #result["FDA_drugs"] =  list(set(IN_list))
    #result["not_FDR_drugs"] = list(set(OUT_list))
    #list("FDA_drugs" = unique(IN_list),"not_FDR_drugs" = unique(OUT_list))
    return(result)

def plot_drug_FDA_annotation(annotation_col_1, colors,data_matrix_annotation):
    states = sorted(list(set(annotation_col_1['States'])))
    dic_state_samples = {}
    dic_state_ratio = {}
    ylist = list()
    objects = states
    for state in states:
        dic_state_samples[state] = FDA_drugs_list(data_matrix_annotation.loc[annotation_col_1[annotation_col_1['States'] == state].index.values]['pert_iname'].values)
        dic_state_ratio[state] = dic_state_samples[state][dic_state_samples[state]['FDA_approved'] == 'In'].shape[0]/dic_state_samples[state][dic_state_samples[state]['FDA_approved'] == 'Out'].shape[0]
        ylist.append(dic_state_ratio[state])
    return(ylist)

def exclude_list(list1, list2):
    result_list = []
    for item in list1:
        if item not in list2:
            result_list.append(item)
    return(result_list)