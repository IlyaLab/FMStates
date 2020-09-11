FROM jupyter/datascience-notebook:latest

RUN conda install -c conda-forge --quiet --yes \
    'r-biocmanager=1.30.10' \
    && \
    R -e 'BiocManager::install(version = "3.10")' && \
    R -e 'BiocManager::install(c("ConsensusClusterPlus", "pheatmap", "bezier", "GSVA", "GSEABase"))' && \
    pip install GSVA==1.0.6 && \
    pip install matplotlib-venn && \
    pip install google-cloud && \
    pip install google-cloud-bigquery && \
    pip install lifelines && \
    pip install plotly==4.5.0 && \
    pip install cyjupyter && \
    pip install torch && \
    pip install pytorch-lightning

WORKDIR /project
CMD [ "/bin/bash" ]
