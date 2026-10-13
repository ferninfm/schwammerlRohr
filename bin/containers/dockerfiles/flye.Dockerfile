FROM continuumio/miniconda3
RUN conda create -y -n flye_env --override-channels -c conda-forge -c bioconda python=3.12 flye && conda clean -afy
ENV PATH="/opt/conda/envs/flye_env/bin:${PATH}"
