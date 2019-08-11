# jrflab modules
[![Build Status](https://travis-ci.org/cBioPortal/cbioportal.svg?branch=master)](https://travis-ci.org/jrflab/modules)

## Introduction
This is the implementation of the jrflab pipeline.

## Installation
The easiest way to download this pipeline is to clone the repository.

```
git clone https://github.com/jrflab/modules.git
```

## Dependencies
- An instance of [anaconda](https://www.anaconda.com) or [miniconda](https://conda.io/en/latest/miniconda.html)
- IMB's Platform Load Sharing Facility (LSF) or Oracle's Sun Grid Engine (SGE) for resource management

### Following R Packages
- [see full dependencies](https://github.com/ndbrown6/modules/blob/master/conda_env/jrflab_modules_env.txt)

## Best practices
	
### Conventions
- BAM files are stored in bam/
- Patient and sample names are coded in the `samples.yaml`
- Sample names cannot contain "_"
- FASTQ files are coded in the `sample.fastq.yaml`
- FASTQ files can be stored in a top-level directory of any name as long as the full path is given in the `sample.fastq.yaml`
- Project specific configurations are stored in the `project_config.yaml`

### Whole genome, whole exome and targeted sequencing
- Alignment:
	* Broad's GATK best practices workflow for handling sequence data
- Mutation calling:
	* Haplotype Caller
	* Platypus
	* MuTect
	* Strelka
	* Scalpel
	* Lancet
- Copy number:
	* FACETS
	* ASCAT
- Annotation:
	* Annovar
	* snpEff
	* SIFT
	* pph2
	* vcf2maf
	* VEP
	* OncoKB
	* ClinVar
- HLA Typing
	* HLA Polysolver

### RNA transcriptome sequencing
- FPKM/RPKM:
	* Tophat
	* Cufflinks
	* STAR
- Fusion detection:
	* Tophat-fusion
	* deFuse
	* Fusion-catcher
- Annotation:
	* OncoFuse

## Detailed usage
[wiki](https://github.com/jrflab/modules/wiki)

## Known issues
Under development

### Known bugs
Under development
