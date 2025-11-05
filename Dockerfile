# Use the Ensembl VEP image
#FROM ensemblorg/ensembl-vep:latest AS vep

# Use a Miniconda base image
FROM continuumio/miniconda3

# Set working directory in the container
WORKDIR /nf-postgwas

# Copy only environment files first (better layer caching)
COPY postgwas_py39.yaml /nf-postgwas/postgwas_py39.yaml
COPY postgwas_py27.yaml /nf-postgwas/postgwas_py27.yaml

# Install necessary packages for other installations (add R + build deps)
ENV DEBIAN_FRONTEND=noninteractive

# Make apt robust over HTTPS (small layer)
RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-transport-https ca-certificates software-properties-common \
    && rm -rf /var/lib/apt/lists/*

# Main system deps (incl. R + graphics/font dev libs for tidyverse/ragg/systemfonts)
RUN apt-get update && apt-get install -y --no-install-recommends \
    git wget unzip cmake g++ cpanminus curl make \
    libmariadb-dev-compat zlib1g-dev libbz2-dev libssl-dev libcurl4-openssl-dev \
    libncurses5-dev libncursesw5-dev libdbi-perl libdbd-mariadb-perl \
    libxml-libxml-perl libxml-parser-perl liblzma-dev \
    r-base r-base-dev gfortran libgmp-dev libmpfr-dev libxml2-dev \
    libblas-dev liblapack-dev xz-utils \
    pkg-config \
    libcairo2-dev \
    libfreetype6-dev libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev \
    libpng-dev libjpeg-dev libtiff5-dev \
    fonts-dejavu-core \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y graphviz && rm -rf /var/lib/apt/lists/*

# --- MAGMA: install binary into image (keep data files in external databases mount) ---
RUN apt-get update && apt-get install -y --no-install-recommends libgomp1 && rm -rf /var/lib/apt/lists/*
ARG MAGMA_URL="https://vu.data.surfsara.nl/index.php/s/zkKbNeNOZAhFXZB/download"
RUN mkdir -p /opt/magma && \
    curl -L -o /tmp/magma_v1.10.zip "${MAGMA_URL}" && \
    unzip /tmp/magma_v1.10.zip -d /opt/magma && \
    if [ -f /opt/magma/magma ]; then chmod 755 /opt/magma/magma && install -m 0755 /opt/magma/magma /usr/local/bin/magma; else echo "MAGMA binary not found in archive"; ls -la /opt/magma; exit 1; fi && \
    rm -f /tmp/magma_v1.10.zip

# ----------------------------
# Global R: install packages
# ----------------------------
# remotes first (use HTTPS CRAN)
RUN R -q -e "install.packages('remotes', repos='https://cloud.r-project.org')"

# CRAN packages (non-Bioconductor) — use HTTPS CRAN
RUN R -q -e "install.packages(c( \
    'tidyverse','magrittr','data.table','patchwork','ggrepel','ggplot2','corrplot','reshape2','dplyr', \
    'haploR','stringr','reticulate','susieR','parallel','R.utils','glmnet','tidyr', \
    'foreach','coloc','gwasrapidd','enrichR','gprofiler2','circlize','heatmaply','RColorBrewer', \
    'easyGgplot2','ggpubr','splitstackshape','GeneOverlap','fuzzyjoin','purrr','topr','optparse', \
    'igraph', 'fields', 'here'), \
    repos='https://cloud.r-project.org')"

# Bioconductor packages
RUN R -q -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); \
    BiocManager::install(c('GenomicRanges','biomaRt','ComplexHeatmap'), ask=FALSE, update=FALSE)"

# GitHub packages
RUN R -q -e "remotes::install_github('oyhel/vautils')"
RUN R -q -e "remotes::install_github('PhanstielLab/Sushi')"
RUN R -q -e "remotes::install_github('annahutch/corrcoverage')"
#RUN R -q -e "remotes::install_github('ZikunY/CARMA')"
#RUN R -q -e "remotes::install_github('PheWAS/PheWAS')"

# Specific RcppEigen version
RUN R -q -e "remotes::install_version('RcppEigen', version='0.3.3.9.3', repos='https://cran.r-project.org')"
RUN R -q -e "remotes::install_version('dbplyr', version='2.2.1', upgrade='never')"

# hyprcoloc (no vignettes)
#RUN R -q -e "options(ask=FALSE); remotes::install_github('jrs95/hyprcoloc', build_vignettes=FALSE, upgrade='never')"

# Installing Java for Nextflow (nf-test requires Java too)
RUN wget https://download.oracle.com/java/21/latest/jdk-21_linux-x64_bin.deb && \
    apt-get install -y ./jdk-21_linux-x64_bin.deb && \
    rm jdk-21_linux-x64_bin.deb

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv /nf-postgwas/nextflow /usr/local/bin/

# ----- Install nf-test (CLI into /usr/local/bin) -----
RUN wget -qO- https://get.nf-test.com | bash && \
    install -m 0755 /nf-postgwas/nf-test /usr/local/bin/nf-test

# Install required Perl modules for VEP
#RUN cpanm Test::Warnings

# Install required Perl modules for VEP and additional ones
# RUN cpanm DBI Archive::Zip \
# Set::IntervalTree \
# JSON \
# PerlIO::gzip \
# Bio::DB::BigFile \
# DBD::mysql || \
# Test::Warnings || \
# cpanm --force IPC::Run Bio::Root::Version XML::DOM

# Install VEP
# RUN mkdir -p /opt/vep && \
# cd /opt/vep && \
# git clone https://github.com/Ensembl/ensembl-vep.git && \
# cd ensembl-vep && \
# perl INSTALL.pl
# perl INSTALL.pl -a cf -s homo_sapiens -y GRCh38 --NO_UPDATE

# External tools (cloned at build time)
#RUN git clone https://github.com/bulik/ldsc.git /nf-postgwas/nf-nf-postgwas/software/ldsc
#RUN git clone https://github.com/FinucaneLab/pops.git /nf-postgwas/nf-nf-postgwas/software/pops
#RUN git clone https://github.com/perslab/depict.git /nf-postgwas/nf-nf-postgwas/software/depict
#RUN git clone https://github.com/JonJala/mtag /nf-postgwas/nf-nf-postgwas/software/mtag
#RUN git clone https://github.com/getian107/PRScs.git /nf-postgwas/nf-nf-postgwas/Software/PRScs


# Install Plink 1.9
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip && \
    unzip plink_linux_x86_64_20210606.zip && \
    mv plink /usr/local/bin/ && \
    rm plink_linux_x86_64_20210606.zip

# Install Plink 2
RUN wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20241124.zip && \
    unzip plink2_linux_x86_64_20241124.zip && \
    mv plink2 /usr/local/bin/ && \
    rm plink2_linux_x86_64_20241124.zip

COPY . /nf-postgwas

# Python 2.7 environment specifically for DEPICT, MTAG and LDSC
RUN conda env create -f /nf-postgwas/postgwas_py27.yaml

# Python 3.9 environment for the rest of the python scripts
RUN conda env create -f /nf-postgwas/postgwas_py39.yaml

# Auto-activate postgwas_py39 for every interactive shell
RUN echo 'source /opt/conda/etc/profile.d/conda.sh && conda activate postgwas_py39' >> /etc/bash.bashrc

# Run the container using bash; make Bash a login shell so 'conda activate' works
SHELL ["/bin/bash", "--login", "-c"]

CMD ["bash", "-l", "-c", "conda activate postgwas_py39 && exec bash"]
