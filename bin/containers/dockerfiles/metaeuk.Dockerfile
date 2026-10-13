FROM continuumio/miniconda3
RUN conda install -y --override-channels -c conda-forge -c bioconda metaeuk && conda clean -afy
ENV PATH="/opt/conda/bin:${PATH}"
