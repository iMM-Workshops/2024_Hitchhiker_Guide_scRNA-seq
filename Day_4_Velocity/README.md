# Velocity

This session is handled by Zhisong He

This session will cover the following topics:
* Trajectory analysis
* Pseudotime analysis
* Fate probability estimation
* RNA velocity analysis

It will include the theoretical session, which will be in the morning, and the practical session, which will be largely in the afternoon. For the practical session, it will needs both R and Python. For that, it is recommended to use `miniconda` or `anaconda` to set up a conda environment. For Windows user, it is recommended to install a WSL (Windows Subsystem for Linux), as it has much better conda support than Windows. For MacOS user, please make sure you have set up the compilation toolkit as mentioned in this [link](https://mac.r-project.org/tools/).</br>

## Table of Content
* [Set up the conda environment](#set-up-the-conda-environment)
* [Example dataset](#example-dataset)
* [Markdown vignette](#markdown-vignette)
* [Slides](#slides)
* [More vignettes](#more-vignettes)

## Set up the conda environment
One suggested way to set up the environment is locally via `conda` to create one conda environment including both Python and R, as well as the packages in need, and then use Jupyter Lab as IDE to run the codes.

Here is a guideline to set up the conda environment. At terminal, do the following:
```
conda create -n env_hitchhiker2024 python=3.9 r-base=4 jupyterlab r-reticulate r-irkernel r-devtools scanpy scvelo cellrank python-igraph r-Seurat=5 cython gsl udunits2 -c conda-forge --solver=libmamba
conda activate env_hitchhiker2024
conda install -c bioconda -c conda-forge velocyto.py r-anndata
```

If you are a Linux/WSL user, then do
```
conda install -c conda-forge gcc gxx
conda install -c bioconda -c conda-forge bioconductor-destiny
```
If you are a MacOS user, or the second command failed, do the following
```
conda install -c conda-forge r-biocmanager cmake
echo 'BiocManager::install("destiny")' | R --vanilla
```

Finally, do the following to install the last R packages
```
echo 'devtools::install_github("farrellja/URD")' | R --vanilla
echo 'devtools::install_github("mojaveazure/seurat-disk")' | R --vanilla
```

In the repository there is also the `environment_hitchhiker2024.yml` file which can be used to create the conda environment. However, it doesn't cover the R packages which are not installed directly with `conda`, and it might fail if you are not using Linux.

For this session, Jupyter Lab is a preferred IDE for this practical session, as one can switch between R and Python kernels easily in the notebook. However, feel free to use any other options (e.g. Rstudio for R, PyCharm for Python, VSCode).

## Standard R configuration
Alternatively, you can just use an R installed with standard way, potentially plus your preferred IDE (e.g. Rstudio). For that, you can do the following at your R session:
```
install.packages('Seurat')
install.packages('BiocManager')
BiocManager::install('destiny')
install.packages('devtools')
devtools::install_github("farrellja/URD")
devtools::install_github("mojaveazure/seurat-disk")
```

## Google Colab
For the Python part, you can also use Google Colab (https://colab.research.google.com/). Once a new notebook is started, you can install the required Python packages with `pip`:
```
!pip install scanpy scvelo cellrank python-igraph leidenalg
```
You can also download a file to the virtual machine so that you can load it into the session. For example:
```
!wget -O 'DS1.h5ad' 'https://polybox.ethz.ch/index.php/s/bUYZE6qPgROBggH/download'
```

## Example dataset
The example dataset used in the practical session is one 10x scRNA-seq experiment included in the following publication: [Organoid single-cell genomic atlas uncovers human-specific features of brain development](https://www.nature.com/articles/s41586-019-1654-9). The Seurat object can be downloaded via this [link](https://polybox.ethz.ch/index.php/s/bjNnfD9I3rwpjIt). It includes the processed RNA assay, plus the "spliced" and "unspliced" assays with the spliced and unspliced transcript count matrices, which are needed for the RNA velocity analysis.

The example dataset is also available as an h5ad file (via this [link](https://polybox.ethz.ch/index.php/s/bUYZE6qPgROBggH)) which can be loaded directly to Python. This is to make sure that people failed to convert the Seurat object to AnnData in R due to the environment setup would still be able to continue in the Python session. It contains only the shared features in the standard Cell Ranger output, the spliced transcript count matrix as well as the unspliced count matrix. The X matrix is the log-normalized count matrix. In its `layers` slot, there are `counts`, `spliced` and `unspliced` count matrices, with the latter two critical for running scVelo.

## Markdown vignette
In this repository there is a Markdown tutorial which includes the codes for the practical session: [pseudotime_trajectory_velocity.md](https://github.com/iMM-Workshops/2024_Hitchhiker_Guide_scRNA-seq/blob/main/Day_4_Velocity/pseudotime_trajectory_velocity.md).

## Slides
The slides of the practical session have been uploaded ([pseudotime_velocity_practical.pdf](https://github.com/iMM-Workshops/2024_Hitchhiker_Guide_scRNA-seq/blob/main/Day_4_Velocity/slides/pseudotime_velocity_practical.pdf)). The slides of the theoretical session is also accessible ([pseudotime_velocity_theory.pdf](https://github.com/iMM-Workshops/2024_Hitchhiker_Guide_scRNA-seq/blob/main/Day_4_Velocity/slides/pseudotime_velocity_theory.pdf)) It includes both the environment setup guideline and the codes. Please note that modifications may be made in the final slides.

## More vignettes:
There are more vignettes available in Github for scRNA-seq data analysis:
* [Tutorial of scRNA-seq data analysis in R](https://github.com/quadbio/scRNAseq_analysis_vignette/blob/master/Tutorial.md)
* [Tutorial of condition comparison analysis of scRNA-seq data](https://github.com/quadbio/scRNAseq_comparison_vignette/blob/master/Tutorial.md)
* [Tutorial of single-cell RNA-ATAC multiomic sequencing data analysis in R](https://github.com/quadbio/scMultiome_analysis_vignette/blob/main/Tutorial.md)
