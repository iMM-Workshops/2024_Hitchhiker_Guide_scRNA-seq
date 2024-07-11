# Adaptive Immune Receptor (AIR)

This session is handled by Lisa Dratva.

If you have any questions or comments regarding the session or future VDJ analysis, please don't hesitate to reach me at [LMD76@cam.ac.uk](mailto:lmd76@cam.ac.uk) and connect on [LinkedIn](https://www.linkedin.com/in/lisadratva/)!

## Installation

0. Python installation

You need to have Python 3.8 or newer installed on your system. If you don't have
Python installed, I recommend installing [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

1. Create a new conda environment

```
conda create --name imm_air_env python=3.10
conda activate imm_air_env
conda install anaconda::git
```

2. Install Cell2TCR from Github and check out the branch ```db_extension``` to use the latest features such as integrated IEDB.org querying. 

```bash
git clone https://github.com/Teichlab/cell2tcr.git
cd cell2tcr
git fetch origin
git checkout db_extension
pip install . 
# make sure the trailing dot is included in the line above
```

3. Install further useful packages.  

```bash
conda install seaborn 
pip install scirpy
conda install -c bioconda bbknn
```

4. Install [JupyterLab](https://jupyter.org/install) to work in notebooks and add a kernel to your environment

``` bash
conda install -c conda-forge jupyterlab
python -m ipykernel install --user --name imm_air_env
```

Now, let's launch `jupyter lab`. First, `cd` into the directory where you have saved the `IMM_AIR.ipynb` and the datasets to, then run:
``` bash
jupyter lab
```
Select the `imm_air_env` kernel as your environment and you're ready to go. 

## Google Colab
If you are having trouble with installing this on your machine, you can use a Google Colab notebook instead. You will need a Google account for this. 
1. Add the file `IMM_AIR.ipynb` to your Google Drive and double-click on it, which should open it directly in Google Colab.
2. Add the files from Google Drive folder (https://drive.google.com/drive/folders/1Uk6pmMRzpwnjfZobabHChDNMoGAd-aHE?usp=drive_link) to your own Google Drive to access them via Google Colab.
3. Execute the first lines of the Google Colab notebook related to mounting your Google Drive and installing the packages.