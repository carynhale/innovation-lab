# modules
[![Build Status](https://travis-ci.org/cBioPortal/cbioportal.svg?branch=master)](https://travis-ci.org/jrflab/modules)

## Introduction
This is an updated fork of the jrflab [modules](https://github.com/jrflab/modules).

## Installation
The easiest way to download this pipeline is to clone the repository.

```
git clone https://github.com/ndbrown6/modules.git
```

## Dependencies
- An instance of [anaconda](https://www.anaconda.com) or [miniconda](https://conda.io/en/latest/miniconda.html)
- IMB's Platform Load Sharing Facility (LSF) for resource management
- [Full dependencies](https://github.com/ndbrown6/modules/blob/master/conda_env/jrflab_modules_env.txt)

## Best practices
	
### Conventions
- BAM files are stored in `bam/`
- Patient and sample names are coded in the `samples.yaml`
- Sample names cannot contain "_"
- FASTQ files are coded in the `sample.fastq.yaml`
- FASTQ files can be stored in a top-level directory of any name as long as the full path is given in the `sample.fastq.yaml`
- Project specific configurations are stored in the `project_config.yaml`

### Whole genome, whole exome and targeted sequencing
- Alignment:
	* Broad's GATK best practices workflow for handling sequence data ([here](https://software.broadinstitute.org/gatk/best-practices/))
	* Marianas & Watlz for MSK-ACCESS ([here](https://github.com/juberpatel))
	* General purpose workflow to collapse UMI and call consensus duplex reads using fgbio ([here](https://github.com/fulcrumgenomics/fgbio))
- Mutation calling:
	* MuTect v1
	* Platypus
	* Strelka
	* Scalpel
	* Lancet
	* Varscan
- Copy number:
	* FACETS
	* ASCAT
	* CNVkit
- Annotation:
	* Annovar
	* snpEff
	* SIFT
	* pph2
	* vcf2maf ([here](https://github.com/mskcc/vcf2maf))
	* VEP
	* OncoKB ([here](https://github.com/oncokb/oncokb-annotator))
	* ClinVar
- Miscellaneous:
	* HLA Polysolver
	* MSIsensor

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
[wiki](https://github.com/ndbrown6/modules/wiki)

## Known issues
Under development

## See also
- [MSKCC](https://github.com/mskcc)
- [cBioPortal](https://github.com/cBioPortal)
- [OncoKB](https://github.com/oncokb)
- [MSK-ACCESS](https://github.com/msk-access)
- [GRAIL](https://github.com/grailbio)

## Contributors
 - See full list [here](https://github.com/ndbrown6/modules/graphs/contributors)


