# Adaptive Immune Receptor (AIR)

This session is handled by Lisa Dratva

## Installation

0. Python installation

You need to have Python 3.8 or newer installed on your system. If you don't have
Python installed, I recommend installing [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

1. Create a new conda environment

```
conda create --name imm_air_env python=3.10
conda activate imm_air_env
```

2. Install Cell2TCR from Github and check out the branch ```db_extension``` to use the latest features such as integrated IEDB.org querying. 

```bash
git clone https://github.com/Teichlab/cell2tcr.git
cd cell2tcr
git fetch origin
git checkout db_extension
pip install .
```

3. Install [JupyterLab](https://jupyter.org/install) to work in notebooks

```
mamba install -c conda-forge jupyterlab
# python -m ipykernel install --user --name imm_air_env # TBD
```

