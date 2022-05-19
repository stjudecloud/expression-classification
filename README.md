# Interactive t-SNE

Code used to generate interactive t-SNE plots using RNA-Seq counts data.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

To get started, you'll want to have Python 3 and R 3.6.1.

```bash
conda create -n itsne-dev \
             -c conda-forge \
             -c bioconda \
             r-base=3.6.1 \
             python=3.7 \
             logzero \
             numpy \
             pandas \
             scikit-learn \
             bioconductor-biocinstaller \
             bioconductor-sva \
             bioconductor-deseq2 \
             r-tsne r-getopt r-plotly r-optparse r-rtsne r-pracma r-data.table -y
conda activate itsne-dev
```

Alternatively, you can install the anaconda dependencies directly from the `environment.yml` file.

```bash
conda env create -f environment.yml
```

### Installing

First, the author recommends you run the following command line script to ensure all R
packages are loaded before your first run:

```bash
Rscript scripts/itsne-normalize-matrix.R
```

From here, you can install the python package:

```bash
python3 setup.py install
```

Download the relevant gene model:

```bash
curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz -o gencode.v31.annotation.gtf.gz
```

And try out the script:

```bash
itsne-main counts/* \
           -b reference/gene.excludelist.tsv \
           -c reference/covariates.tsv \
           -g gencode.v31.annotation.gtf.gz \
           -o tsne.html
```

## Running the tests

No tests currently exist. If and when tests are added, the authors will fill in this section.

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/stjudecloud/expression-classification/tags). 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

