# Data Science Capstone Fall 2022
## Authors: Ju, Katharina, Darvesh

# Setting Up
Most files that will be needed are contained within the repository. The only files that are not contained within the repository are those files in `raw_data/Pasca_scRNAseq/inter/`. They are not included because they are too large to commit to GitHub via normal means and we lack adequate storage space and bandwidth on git lfs (Git Large-file Storage). However, these files are also generated via running the pipeline so they are not strictly necessary for running the code from scratch. A `git clone` should be sufficient to get started. If you have access to the Google Drive folder with the data that belongs in `raw_data/Pasca_scRNAseq/inter/`, go ahead and copy it there. The files in there will automatically be ignored via the `.gitignore` file so any changes to them will not be tracked. This may change in the future but it is unlikely due to the extremely large nature of the files (>16 Gb).

## Python
This repo uses python 3.10.4 for all the python code (scripts and notebooks). We recommend creating a conda environment as follows:

```
conda env create -f environment.yml
```

This will install all the necessary python dependencies on your machine. To access this environment, type the following (`scrna` is the default name of the conda environment):

```
conda activate scrna
```

## R
We use version 4.2.1 for R. The R packages are installed in the scripts so there's no need to worry about installing the R packages. This may change later. Please download the appropriate version of R [here](https://archive.linux.duke.edu/cran/) based on your operating system. 

## Docker
Coming soon...