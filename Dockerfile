FROM ubuntu:18.04 as builder
RUN apt-get update \
    && apt-get upgrade -y \ 
    && apt-get --yes install \
        build-essential \
        openjdk-11-jdk-headless \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" -O miniconda.sh && \ 
    /bin/bash miniconda.sh -b -p /opt/conda/ && \
    rm miniconda.sh

ENV PATH /opt/conda/bin:$PATH

RUN conda update -n base -c defaults conda -y && \
    conda install \
    -c conda-forge \
    -c bioconda \
             r-base=3.6.1 \
             python=3.8 \
             pandas=1.1.0 \
             bioconductor-biocinstaller \
             bioconductor-sva \
             bioconductor-deseq2 \
             scikit-learn \
             r-tsne r-getopt r-plotly  \
             pandoc -y \
    && conda clean --all -y

RUN Rscript -e 'install.packages(c("optparse", "Rtsne", "plotly", "sva", "stringr", "pracma"), dependencies=TRUE, repos="http://cran.rstudio.com/")'

RUN Rscript -e 'install.packages("BiocManager", dependencies=TRUE, repos="http://cran.rstudio.com/")' 
RUN Rscript -e 'BiocManager::install(version = "3.10")'
RUN Rscript -e 'BiocManager::install("DESeq2", update = TRUE, ask = FALSE)'

COPY scripts /opt/tsne/scripts
COPY itsne /opt/tsne/itsne
COPY setup.py /opt/tsne
COPY requirements.txt /opt/tsne
COPY environment.yml /opt/tsne

WORKDIR /opt/tsne
RUN  python3 setup.py install

WORKDIR /
