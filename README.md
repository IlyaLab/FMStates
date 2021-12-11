Usage: 
Before building the framework, make sure docker is installed and running in your system. All datasets used in the manuscript of "Functional module states framework reveals cell states for drug and target prediction" by Guangrong Qin et al., can be found in https://osf.io/34xnm/?view_only=5b968aebebe14d4c97ff9d7ce4cb5070

After clone all the scripts in the repo, the following steps will help you to get familar with the notebooks developed in this project. 

Step 1: imitates and build the environment
./cli.rb build

Step 2: run the docker image
./cli.rb run-notebook

It will initiate the notebook server in a docker environment, to access the notebook, open the URLs in a browser. 

Example jupyter notebooks:
To run the example jupyter notebook as shown in the manuscript, download all files from the directory of 'Sample_input' from https://osf.io/34xnm/?view_only=5b968aebebe14d4c97ff9d7ce4cb5070, and put it under the 'project' directory.

Example 1: Functional transcriptional states corresponding to drug treatment.
### Generate module specific Transcription factors (TFs) with enriched target genes in the context
```
$ Example1-generate-TF-pairs.ipynb
```
### Generate functional module matrix
```
$ Example1-generate-FM-matrix.ipynb
```
###  Compare FM factor matrix between two groups to define the relative FM-factors, and use them for clustering (define states) and annotation.
```
$ Example1-comparing-clustering-annotation.ipynb
```
### Annotate each stats using module specific TFs

```
$ Example1_annotation_tf.ipynb 
```
Users can use Cytoscape for visualizing the regulation of the transcription factors and the regulated pathways.

### Annotate each state using external data: drug response (whether different drug dosages are associated with different states)
```
$ Example1_annotation_drugResponse.ipynb
```

### Annotate each state using external data: drug targets (whether the drug target classes are associated with different states)
```
Example1_annotation_targets.ipynb
```


Example 2: Predict targetable vulnerability using the transcriptomic states

### Generate FM-factors for the transcriptome of cell lines after shRNA knockdown.
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
Example3_association_analysis_FMfactor_drugResponse.ipynb
```

### Predict drug response using random forest.
```
Example3_predict_drug_response_rf.ipynb
```

Example 4: Applying the FM-States framework to patient derived samples in AML
### FM-factors generation, clustering, association analysis and prediction models are imbedded in the following script
Example4_AML_analysis.ipynb
###

