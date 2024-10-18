
# PyViscount

A Python framework for validation of false discovery rate (FDR) estimation methods in shotgun proteomics using random search space partition.


## Features

- Two validation modes: pre- and post-search
- Compatible with output of 4 popular proteomics search engines
- Two reporting modes: contour and optimized single-line plots


## Installation

It is recommended to use PyViscount in a conda environment. To create and activate a new environment with all necessary Python packages (provided in "requirements.txt"), run the following commands:

```bash
conda create --name <env> --file requirements.txt
conda activate <env>
```

Then, clone the PyViscount repository to your local machine:

```bash
git clone https://github.com/dommad/pyviscount.git
cd pyviscount
```

Finally, you can install the package:

```bash

pip install .
```

Alternatively, you can install the package in editable mode, which allows you to make changes to the source code without needing to reinstall the package:

```bash
pip install -e .
```

To verify that PyViscount was installed successfully, you can run:

```bash
pip show pyviscount
```

or within a Python environment:

```python
import pyviscount
```
## Input data

PyViscount accepts the output provided by 4 search engines: Tide, Comet, MSFragger, and SpectraST. Allowed file formats: .txt, .pep.xml, .pepXML, .mzid.

PyViscount requires a configuration file (example [here](https://github.com/dommad/pyviscount/blob/main/data/config.ini)). Some of the important options to be specified in the configuration file:

- *mode* (mode of analysis; postsearch or presearch)
- *engine* (engine used to generate the input files)
- *fdr_score* (FDR is calculated with respect to this score)
- *threshold_score* (quality filtering is conducted using this score as a cutoff)
- *output_path* (path to directory in which all output files are placed)


## Usage

PyViscount can be run via CLI using the following command:

```bash

pyviscount -conf configuration.ini -t target_file.txt -td target_decoy_file.txt -d decoy_file.txt

```

It can also be used in interactive Python environments (like Jupyter notebook):

```python

from pyviscount import run

config_file = "./config.ini"
target = './tide-search.target.txt'
target_decoy = './tide-search.txt'
decoy = './tide-search.decoy.txt'

results = run.run_postsearch_validation(config_file, target, target_decoy, decoy)

```

You can find a sample Jupyter notebook [here](https://github.com/dommad/pyviscount/blob/main/examples/test.ipynb).

