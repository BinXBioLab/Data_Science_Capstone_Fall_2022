FROM ubuntu:22.10

ENV DEBIAN_FRONTEND noninteractive
ENV CONDA_DIR /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH 

RUN apt-get -y update
RUN apt-get install -y python3.10 wget git

# Installing R 
RUN apt update -qq
RUN apt install --no-install-recommends -y software-properties-common dirmngr
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/'
RUN apt install --no-install-recommends -y r-base

# Installing anaconda and python dependencies
# Install miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && /bin/bash ~/miniconda.sh -b -p /opt/conda
RUN conda update -n base -c defaults conda
RUN git clone https://github.com/BinXBioLab/Data_Science_Capstone_Fall_2022.git
RUN cd Data_Science_Capstone_Fall_2022
RUN pip install -r requirements.txt

