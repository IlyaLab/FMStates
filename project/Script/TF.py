import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
import scipy 
import FM_States
import scipy.stats as ss
import statsmodels
from statsmodels import stats
from statsmodels.stats import multitest


#### User input: User can select functional moduels and input the gene expression matrix 

ROOT_DIR = os.path.abspath("../")

### Transcription factors
#### Step 1. select transcription factor and target gene pairs with different evidences
def sele_pairs(evidence):
    ## evidence should be a list, among the list, only 'chipSeq','TFbindingMotif','curated' are effective selection.
    TF_pair_all = pd.read_csv(os.path.join(ROOT_DIR, 'Dataset/database.csv')) #Data from Luz Garcia-Alonso et al., Genome biology, 2019. 
    #TF_pair_all = pd.read_csv("https://genome.cshlp.org/content/suppl/2021/03/02/gr.240663.118.DC2/Revised_Supplemental_Table_S3.csv")
    if 'chipSeq' in evidence:
        TF_pair_chipSeq = TF_pair_all.loc[TF_pair_all['is_evidence_chipSeq'] == True]
    else:
        TF_pair_chipSeq = pd.DataFrame()

    if 'TFbindingMotif' in evidence: 
        TF_pair_TFbindingMotif = TF_pair_all.loc[TF_pair_all['is_evidence_TFbindingMotif'] == True]
    else:
        TF_pair_chipSeq = pd.DataFrame()

    if 'curated' in evidence:
        TF_pair_curated = TF_pair_all.loc[TF_pair_all['is_evidence_curateddatabase'] == True]
    else:
        TF_pair_curated = {}

    result = pd.concat([TF_pair_curated, TF_pair_chipSeq,TF_pair_TFbindingMotif])
    return(result)

#### 2. select transcription factor and target gene pairs with co-expression evidences
def pair_cor_test(test, data_matrix_MCF7_CTRP2):
    from scipy.stats import pearsonr
    TF_remain = []
    Target_remain = []
    corr_list = []
    dic_TF_target_pair = {}
    for i in range(0,test.shape[0]):
        dic_TF_target_pair[test.iloc[i]['TF'] + ':' + test.iloc[i]['target']] = '' 
    #test = x_A.loc[x_A['target'].isin(dic_module[list(set(KEGG_level2['name']))[i]])]
    
    for TF in list(set(test['TF'])):
        for Target in list(set(test['target'])):
            if (TF+':'+Target) in dic_TF_target_pair and (TF in data_matrix_MCF7_CTRP2.columns.tolist()) and (Target in data_matrix_MCF7_CTRP2.columns.tolist()) :
                corr, p = pearsonr(data_matrix_MCF7_CTRP2[TF],data_matrix_MCF7_CTRP2[Target])
                
                if p < 0.05 and (abs(corr)) > 0.2:  ##abs(corr) > 0.2 for both activator and supressor or corr > 0.2 for only activator
                
                    TF_remain.append(TF)
                    Target_remain.append(Target)
                    corr_list.append(corr)

    pair_cor = pd.DataFrame({"TF":TF_remain,"Target":Target_remain,"corr":corr_list})
    return(pair_cor)

#### 3. select transcription factor and target gene pairs with both co-expression evidence and high confident evidences from the orignal literatures.
def get_tf_target_pair(data_matrix_MCF7_CTRP2,pathway_select,dic_module):
    pair_cor_df = pd.DataFrame({})
    for i in range(0,len(pathway_select)):
        x = sele_pairs(['chipSeq','TFbindingMotif','curated']) 
        x_AB = x.loc[x['score'].isin(['A','B'])]
        x_AB_module = x_AB.loc[x_AB['target'].isin(dic_module[pathway_select[i]])]
        pair_cor = pair_cor_test(x_AB_module, data_matrix_MCF7_CTRP2)
        pair_cor_df = pd.concat([pair_cor_df, pair_cor])
        
    Resulting_pairs = pair_cor_df.drop_duplicates()
    return(Resulting_pairs)

#### 4. Estimate which transcription factors have the key regulation role for specific module.
def get_tfpairs_for_select_pathways(data_matrix_MCF7_CTRP2,pathway_select,dic_module): ## remove parameter active = True
    
    Resulting_pairs = get_tf_target_pair(data_matrix_MCF7_CTRP2,pathway_select,dic_module)
    #if active == True:
    #    Resulting_pairs = Resulting_pairs.loc[Resulting_pairs['corr'] > 0]
    #elif active == False:
    #    Resulting_pairs = Resulting_pairs.loc[Resulting_pairs['corr'] < 0]
        
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


def High_req_TFs(test):
    tem_list_uniq = list(set(test['TF']))
    count_list_uniq = []
    tem_list = test['TF'].tolist()
    for tf in tem_list_uniq:
        tem_len = tem_list.count(tf)
        count_list_uniq.append(tem_len)
    y = pd.DataFrame({'TF':tem_list_uniq, 'Count':count_list_uniq })
    y_sort = y.sort_values(by = 'Count', ascending=False)

    plt.figure(num=None, figsize=(20, 4), dpi=100, facecolor='w', edgecolor='k')
    plt.bar(y_sort['TF'], y_sort['Count'])
    plt.xticks(y_sort['TF'],  rotation='vertical')
    return(y_sort)

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

def TF_annotation(x,TF_pairs,data_matrix_MCF7_CTRP2, annotation_col_1, output_dir, active = True):
    TF_used_df = TF_pairs[TF_pairs['Occupancy_inPathway'] > 0.03]
    TF_used_df = TF_used_df.loc[TF_used_df['Occupancy_inTF'] > 0.2]
    #TF_used_df = TF_pairs_active
    pathway_used = set()
    for item in x.columns.tolist():
        pathway_used.add(item.split('_')[0])
    pathway_used = list(pathway_used)
    TF_used_df = TF_used_df.loc[TF_used_df['Pathway'].isin(pathway_used)]
    TF_used = set(TF_used_df['TFs'])

    #data_matrix_MCF7_CTRP2 = pd.read_csv(os.path.join(ROOT_DIR, "Sample_input/Example1/Sample1_data_MCF7_drugs_CTRP2.csv"))
    #data_matrix_MCF7_CTRP2 = data_matrix_MCF7_CTRP2.set_index('Unnamed: 0')
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

            cohen_list_tmp.append((FM_States.cohen_dist(t1,t2)))
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


def cytoscape_plot(x, F1, F2,Edge, types):
    from cyjupyter import Cytoscape
    import json

    nodes_list = []
    edges_list = []

    if types == 'category':
        ID1_list = []
        ID2_list = []
        for i in range(0, x.shape[0]):
            ID1_list.append( str(x[F1].values[i]) )
            ID2_list.append( str(x[F2].values[i]) )
            
        x['ID1'] = ID1_list
        x['ID2'] = ID2_list

        for i in range(0, x.shape[0]):
            nodes_list.append({'data':{'id':x['ID1'].values[i], 
                                       'color':'yellow', 
                                       'shape':'ellipse'}} )
            
            nodes_list.append({'data':{'id':x['ID2'].values[i], 
                                       'color':'grey', 
                                       'shape':'rectangle'}} )
            
            if x[Edge].values[i] == 'Active' :
                edges_list.append({'data':{'id': 'edge'+str(i),
                                           'source': x['ID1'].values[i],
                                           'target': x['ID2'].values[i], 
                                           'type': 'red'}})
                
            elif x[Edge].values[i]  == 'Supress':
                edges_list.append({'data':{'id': 'edge'+str(i),
                                           'source': x['ID1'].values[i],
                                           'target': x['ID2'].values[i], 
                                           'type': 'black'}})

    my_cy = {'elements':{'nodes':nodes_list,'edges':edges_list}}

    mystyle = [{
            'selector': 'node',
            'style': {
                'background-color': 'data(color)',
                'label': 'data(id)',
                'width': 12,
                'height': 12,
                'shape':'data(shape)',
                'color': 'grey',
                'font-weight': 400,
                'text-halign': 'middle',
                'text-valign': 'bottom',
                'font-size': 6,
                'size':3
            }},
            {
            'selector': 'edge',
            'style': {
                'width': 1,
                'line-color': 'data(type)',
                'target-arrow-color': '#37474F',
                'target-arrow-shape': 'triangle'}
            }]
    cy = Cytoscape(data = my_cy, visual_style = mystyle,  background =('#FFFFFF'))
    return(cy)
