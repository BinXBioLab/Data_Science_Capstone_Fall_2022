# Data Science Capstone Fall 2022
## Authors: Ju, Katharina, Darvesh

# Setting Up
Most files that will be needed are contained within the repository. The only files that are not contained within the repository are those files in `raw_data/Pasca_scRNAseq/inter/`. They are not included because they are too large to commit to GitHub via normal means and we lack adequate storage space and bandwidth on git lfs (Git Large-file Storage). However, these files are also generated via running the pipeline so they are not strictly necessary for running the code from scratch. A `git clone` should be sufficient to get started. If you have access to the Google Drive folder with the data that belongs in `raw_data/Pasca_scRNAseq/inter/`, go ahead and copy it there. The files in there will automatically be ignored via the `.gitignore` file so any changes to them will not be tracked. This may change in the future but it is unlikely due to the extremely large nature of the files (>16 Gb).

## Python
We have tried to resolve as many conflicts between versions of python as possible but there may be some issues. The best way to install the necessary requirements is as follows:

```
pip install -r base-requirements.txt -r cpdb-requirements.txt -r tangram-requirements.txt
```

or 

```
pip install -r requirements.txt
```

```
conda activate scrna
```

## R
We use version 4.2.1 for R. The R packages are installed in the scripts so there's no need to worry about installing the R packages. This may change later. Please download the appropriate version of R [here](https://archive.linux.duke.edu/cran/) based on your operating system. 

## Docker
Coming soon...