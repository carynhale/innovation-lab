# innovation-lab
[![Build Status](https://travis-ci.com/ndbrown6/innovation-lab.svg?token=WQkjmC3gu8Nd4XcQmFkn&branch=master)](https://travis-ci.com/ndbrown6/innovation-lab)

## Introduction
This is a lightweight fork of a previous project [modules](https://github.com/ndbrown6/modules).

## Installation
For outside users, the easiest way to download/ execute this pipeline is to install
locally miniconda with its dependencies and clone this git repository

```
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O $HOME/miniconda.sh
bash $HOME/miniconda.sh -b -p $HOME/.miniconda
export PATH="$HOME/.miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda info -a
conda config --add channels r
conda config --add channels bioconda
mkdir -p $HOME/share $HOME/share/usr $HOME/share/usr/env
export INNOVATION_ENV_VERSION=0.0.1
conda create -q -p $HOME/share/usr/env/innovation-lab-$INNOVATION_ENV_VERSION
source activate $HOME/share/usr/env/innovation-lab-$INNOVATION_ENV_VERSION
git clone https://github.com/ndbrown6/innovation-lab.git
```

## Dependencies
- An instance of [anaconda](https://www.anaconda.com) or [miniconda](https://conda.io/en/latest/miniconda.html)
- IMB's Platform Load Sharing Facility (LSF) for resource management
- [Full dependencies](https://github.com/ndbrown6/innovation-lab/tree/master/conda)

## Best practices
	
### Conventions

### Whole genome, whole exome and targeted sequencing

## Detailed usage

## Known issues

## See also

## Contributors
