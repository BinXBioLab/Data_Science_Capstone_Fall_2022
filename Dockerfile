FROM ubuntu:22.10

ENV DEBIAN_FRONTEND noninteractive
ENV CONDA_DIR /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH 

RUN apt-get -y update
RUN apt-get install -y python3.9 wget git make

# Installing R 
RUN apt update -qq
RUN apt install --no-install-recommends -y software-properties-common dirmngr
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/'
RUN apt install --no-install-recommends -y r-base

# Installing miniconda 
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && /bin/bash ~/miniconda.sh -b -p /opt/conda
RUN conda update -n base -c defaults conda

# Installing python dependencies
COPY requirements.txt /requirements.txt
RUN pip install -r requirements.txt

# Install R dependencies
RUN R -e "install.packages(c('BiocManager'))"
RUN R -e "BiocManager::install(c('Seurat', 'SingleR'))"
RUN R -e "install.packages(c('Rcpp', 'ggplot2', 'liger', 'reshape2', 'plyr', 'MAST', 'scran', 'gam', 'clusterExperiment'))"
