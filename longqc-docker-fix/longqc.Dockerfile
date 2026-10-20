FROM continuumio/miniconda3

RUN apt-get update && apt-get install -y --no-install-recommends \
    git make gcc g++ zlib1g-dev ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/yfukasawa/LongQC.git /opt/LongQC \
    && make -C /opt/LongQC/minimap2-coverage

# Do not modify the base Conda environment: current Miniconda images may
# depend on a newer Python ABI than the one required by LongQC.
RUN conda create -n longqc -y --override-channels \
      -c conda-forge -c bioconda \
      python=3.9 numpy pandas scipy jinja2 h5py matplotlib \
      scikit-learn edlib pysam python-edlib \
    && conda clean -afy

ENV PATH="/opt/LongQC:/opt/conda/envs/longqc/bin:/opt/conda/bin:${PATH}"
