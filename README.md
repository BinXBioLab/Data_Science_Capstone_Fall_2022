# Data Science Capstone Fall 2022
## Authors: Ju, Katharina, Darvesh

# Setting Up
Most files that will be needed are contained within the repository. The only files that are not contained within the repository are those files in `raw_data/Pasca_scRNAseq/inter/`. They are not included because they are too large to commit to GitHub via normal means and we lack adequate storage space and bandwidth on git lfs (Git Large-file Storage). However, these files are also generated via running the pipeline so they are not strictly necessary for running the code from scratch. A `git clone` should be sufficient to get started. If you have access to the Google Drive folder with the data that belongs in `raw_data/Pasca_scRNAseq/inter/`, go ahead and copy it there. The files in there will automatically be ignored via the `.gitignore` file so any changes to them will not be tracked. This may change in the future but it is unlikely due to the extremely large nature of the files (>16 Gb).

## Local Machine - With Docker
### Preprocessing & Main Pipeline 
To run the preprocessing and main pipeline steps with Docker on your local computer, you will first need to build a Docker image. Docker is a containerization tool used to create isolated environments that are portable across operating systems (i.e. Windows, MacOS, and Linux). Dockerfiles define the specification for how to create an image and each image can be run as a container. One image can be used for multiple containers. 

1. Install Docker on your computer ([installation instructions are here](https://docs.docker.com/get-docker/)). 
2. Clone this repository to your local computer via 

```
git clone git@github.com:BinXBioLab/Data_Science_Capstone_Fall_2022.git
```

3. Make sure your working directory is the root of the repository you cloned. If you cloned with the command above you can type the following to change into the root directory of the repository

```
cd ./Data_Science_Capstone_Fall_2022
```

4. Build the Docker image for the preprocessing and pipeline steps via the command below. This step can take a long time specifically because of the `Seurat` and `SingleR` installation.

```
docker build -f Dockerfile-Pipeline -t preprocessing/pipeline .
```

5. To test whether the container was successfully built you can run the command below.

```
docker run --name preprocessing --rm -v $PWD/optimized_codes:/optimized_codes -i -t preprocessing/pipeline python /optimized_codes/scripts/preprocessing.py --help
```

6. Check the output of step 5 to make sure it looks like this. This is the expected output of the `--help` flag on preprocessing.
```
Usage: preprocessing.py [OPTIONS]

  Top-level function which executes all the high-level functions

  Args:     in_dir (str): Root of directory with relevant inputs     inter_dir
  (str): Directory path to store intermediate files     pct_count_lower (int):
  lower bound for percent count filtering     pct_count_upper (int): upper
  bound for percent count filtering     genes_by_count_lower (int): lower
  bound for gene count filtering     genes_by_count_upper (int): upper bound
  for gene count filtering     total_counts_lower (float): lower bound for
  total count filtering     total_counts_upper (float): upper bound for total
  count filtering     aggregated_filename (str): Filename for aggregated
  annotated data object when saving to disk

Options:
  --in_dir TEXT                   Directory to look for 10x Genomics
                                  directories
  --inter_dir TEXT                Directory to store intermediate outputs
  --pct_count_lower TEXT          Lower bound of percent counts to filter
  --pct_count_upper INTEGER       Upper bound of percent counts to filter
  --genes_by_count_lower INTEGER  Lower bound of genes by counts to filter
  --genes_by_count_upper INTEGER  Upper bound of genes by counts to filter
  --total_counts_lower TEXT       Lower bound of total counts to filter
  --total_counts_upper FLOAT      Upper bound of total counts to filter
  --aggregated_filename TEXT      Filename for concatenated counts + labels
                                  h5ad data file
  --help                          Show this message and exit.
```

### Tangram (With Docker)
There is no non-Docker version available at this moment, but you are free to try using Tangram without Docker. The steps to use Tangram with Docker are as follows

1. Install Docker on your computer ([installation instructions are here](https://docs.docker.com/get-docker/)). 

2. Clone this repository to your local computer via 

```
git clone git@github.com:BinXBioLab/Data_Science_Capstone_Fall_2022.git
```

3. Make sure your working directory is the root of the repository you cloned. If you cloned with the command above you can type the following to change into the root directory of the repository

```
cd ./Data_Science_Capstone_Fall_2022
```

1. Build the Docker image for the Tangram with the command below (it's important to be at the root of this repository). This step can take a long time but it shouldn't take longer than the main pipeline + preprocessing Docker container. The container requires approximately **6 Gb** of storage to build.

```Docker
docker build -f Dockerfile-Tangram -t tangram .
```

4. To test whether the container was successfully built you can run the command below.

```Docker
docker run -it --rm --gpus all tangram:latest python -c "import torch; print('Pytorch Version:', torch.__version__); print('NVIDIA GPU Available:', torch.cuda.is_available())"
```

The output should look like this if everything was installed correctly
```python
Pytorch Version: 1.13.0
NVIDIA GPU Available: True
```

If you are on a machine without a GPU but the container built successfully you will see:
```python
Pytorch Version: 1.13.0
NVIDIA GPU Available: False
```

5. To run code utilizing Tangram, you'll first need to create a folder which will be mounted to the container. The purpose of this folder is to share only the files you want between the container and your host machine. This folder will contain things like the input files, figures generated by Tangram, and anything else you want to use in the container during code execution. For example, if you create a folder in your current folder called `data/` which contained `Tangram_CPU.py` you would run the following to run the script within the container.

```Docker
docker run -it --rm --gpus all -v $PWD/data:/home tangram:latest python /home/Tangram_CPU.py
```

`$PWD` is a shorthand for 'print working directory'. This is used because Docker does not allow characters like `.` in its commands. As a consequence something like `./data:/home` would fail to run. If you're specifying a path in a different folder, make sure you write the full path.

You will also need to make sure you account for the mounted destination when reading and writing files in the container. For example, by default python will start in the `/workspace` directory in the container despite the presence of a mounted directory. For example, if you used the `-v $PWD/data:/home` argument in your `docker run` command (like the example above), you could add the following to the top of the script.

```python
import os
os.chdir("/home")
```

Alternatively, you could specify absolute paths for the files like so

```python
...
some_file_path = "/home/path/to/file"
...
```

6. To run the code from the repository use the following command

```
docker run -it --rm --gpus all -v $PWD/data:/home tangram:latest python /home/Tangram_CPU.py
```

