FROM ubuntu:20.04 as builder

RUN apt-get update \
    && apt-get upgrade -y \ 
    && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get --yes install \
        build-essential \
        openjdk-11-jdk-headless \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN wget "https://repo.anaconda.com/miniconda/Miniconda3-py38_23.5.2-0-Linux-x86_64.sh" -O miniconda.sh && \
    /bin/bash miniconda.sh -b -p /opt/conda/ && \
    rm miniconda.sh

ENV PATH /opt/conda/bin:$PATH

RUN conda update -n base -c defaults conda -y && \
    conda install \
    -n base \
    -c conda-forge \
    conda-libmamba-solver && \
    conda config --set solver libmamba

RUN conda install \
    -c conda-forge \
    -c bioconda \
             r-base \
             python=3.8 \
             pandas \
             bioconductor-biocinstaller \
             bioconductor-sva>3.35.2 \
             bioconductor-deseq2 \
             scikit-learn \
             r-rtsne r-getopt r-plotly \
             r-rcppeigen r-hmisc r-optparse \
             r-pracma \
             pandoc -y \
    && conda clean --all -y

ENV R_LIBS_SITE $R_LIBS_SITE:/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library:/opt/conda/lib/R/library/

COPY scripts /opt/tsne/scripts
COPY itsne /opt/tsne/itsne
COPY setup.py /opt/tsne
COPY requirements.txt /opt/tsne
COPY environment.yml /opt/tsne

WORKDIR /opt/tsne
RUN pip install

WORKDIR /

