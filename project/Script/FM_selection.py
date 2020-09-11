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


def Jaccard_coef_simple(list1, list2):
    ''' Compute the Jaccard coeffiency for two list!
        Parameter: 
            list1: 
            list2: 
        Output:
            Jaccard coeffiency. 
    '''
    intersection_num = set(list1).intersection(set(list2))
    union_num = set(list1 + list2)
    coef = len(intersection_num)/union_num
    return(coef)

def Jaccard_coef(module_selected_gmt, name, member):
    namelist = list(set(module_selected_gmt[name]))
    namelist_withNum = []
    result = pd.DataFrame()
    for name1 in namelist:
        num = len(set(module_selected_gmt.loc[module_selected_gmt[name] == name1][member]))
        namelist_withNum.append(name1 + '(' + str(num)+')')
        tem_coef = []
        for name2 in namelist:
            coef = Jaccard_coef_simple( module_selected_gmt.loc[module_selected_gmt[name] == name1][member], module_selected_gmt.loc[module_selected_gmt[name] == name2][member])
            tem_coef.append(coef)
        result[name1] = tem_coef 
    result.index = namelist_withNum
    
    return(result)

### Load functional modules
def load_function_modules(label):
    if label == "KEGG":
        KEGG_level3 = pd.read_csv(os.path.join(ROOT_DIR,"Dataset/KEGG_GENE_level3.csv"),header= None, names = ("name","member"))
        KEGG_level2 = pd.read_csv(os.path.join(ROOT_DIR,"Dataset/KEGG_GENE_level2_select.csv"),header= None, names = ("name","member"))
        KEGG_modules = pd.concat([KEGG_level2, KEGG_level3])

        KEGG_modules['description'] = KEGG_modules['name']
        KEGG_modules.columns = ["name",'member','description']

        dic_module = {}
        module_list = list(set(KEGG_modules['name']))#
        for module in module_list:
            dic_module[module] = list(set(KEGG_modules[KEGG_modules['name'] == module]['member']))
        return(dic_module, KEGG_level2, KEGG_level3, KEGG_modules)
