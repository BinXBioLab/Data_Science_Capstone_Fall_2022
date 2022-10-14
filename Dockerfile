FROM ubuntu:22.10

ENV DEBIAN_FRONTEND noninteractive
ENV CONDA_DIR /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH 

RUN apt update
RUN apt install -y software-properties-common dirmngr
RUN apt update
RUN apt install -y python3 wget git make
RUN apt update

# Installing R 
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/'
RUN apt install --no-install-recommends -y r-base

# Installing python dependencies
COPY requirements.txt /requirements.txt
RUN pip install -r requirements.txt

# Install R dependencies
RUN R -e "install.packages(c('BiocManager'))"
RUN R -e "BiocManager::install(c('Seurat', 'SingleR'))"
RUN R -e "install.packages(c('Rcpp', 'ggplot2', 'liger', 'reshape2', 'plyr', 'MAST', 'scran', 'gam', 'clusterExperiment'))"
