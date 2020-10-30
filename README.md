Usage: 
Before building the framework, make sure docker is installed and running in your system. All datasts used in the manuscript of "Functional module states framework reveals cell states for drug and target prediction" by Guangrong Qin et al., can be found in https://osf.io/34xnm/?view_only=5b968aebebe14d4c97ff9d7ce4cb5070

Step 1: initate and build the environment
./cli.rb build

Step 2: run the docker image
./cli.rb run-notebook

It will initate the notebook server in a docker environment, to access the notebook, open the URLs in a brower. 


Example jupyter notebooks:

Example 1: Functional transcriptional states corresponding to drug treatment.
### Generate content specific Transcription factor (TF) and enriched targeted pathways
```
$ Example1-generate-TF-pairs.ipynb
```
### Generate functional module matrix
```
$ Example1-generate-FM-matrix.ipynb
```
###  Compare FM matrix between two groups to define the relative FM-factors, and use them for clustering (define states) and annotation.
```
$ Example1-comparing-clustering-annotation.ipynb
```
### Annotate each states using TF-regulation network

```
$ Example1_annotation_tf.ipynb 
```
Users can use cytoscape to further visulize the regulation of the transcription factors and the regulated pathways.

### Annotate each states using external data: drug response 
```
$ Example1_annotation_drugResponse.ipynb
```

### Annotate each states uing external data: drug targets (Need revision)
```
Example1_annotation_targets.ipynb
```



Example 2: Predict targetable vulnerability using the transcriptome states

### Generate FM-factors for the transcriptome after sh-knockdown.
```
Example2_generate_FM_matrix.ipynb
 ```

 ### Classify gene knockdown to the states define using the transcriptome after drug treatment.
 ```
Example2_KNN_ALL.ipynb
 ```

 ### Analyze the cancer dependency for genes that induce different states
```
Example2_Depmap-allSet.ipynb
```


Example 3: The association between the transcriptional states prior to drug treatment and drug response 
### Generate FM-factors for breast cancer cell line using the transcriptome before drug treatment. 
```
Example3_generate_FM_matrix.ipynb
```

### Association analysis using the FM-factors with drug response data (IC50)
```
Example3_association_analysis_FM-Facotors_drugResponse.ipynb
```

### Predict drug response using random forest.
```
Example3_predict_drug_response_rf.ipynb
```
