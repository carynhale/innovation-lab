# innovation-lab
[![Build Status](https://travis-ci.com/ndbrown6/innovation-lab.svg?token=WQkjmC3gu8Nd4XcQmFkn&branch=master)](https://travis-ci.com/ndbrown6/innovation-lab)

## Introduction
This is a lightweight fork of the jrflab [modules](https://github.com/jrflab/modules).

## Installation
The easiest way to download this pipeline is to clone the repository.

```
git clone https://github.com/ndbrown6/innovation-lab.git
```

## Dependencies
- An instance of [anaconda](https://www.anaconda.com) or [miniconda](https://conda.io/en/latest/miniconda.html)
- IMB's Platform Load Sharing Facility (LSF) for resource management
- [Full dependencies](https://github.com/ndbrown6/innovation-lab/tree/master/conda)

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
	* [Broad's GATK best practices workflow](https://software.broadinstitute.org/gatk/best-practices/) for handling sequence data
	* [Marianas](https://github.com/juberpatel/Marianas) & [Waltz](https://github.com/juberpatel/Waltz) for [MSK-ACCESS](https://github.com/ndbrown6/modules/wiki/4.-MSK%E2%80%90ACCESS)
	* General purpose workflow to collapse UMI and call consensus duplex reads using [fgbio](https://github.com/fulcrumgenomics/fgbio)

## Detailed usage
[wiki](https://github.com/ndbrown6/innovation-lab/wiki)

## Known issues
Under development

## See also
- [MSKCC](https://github.com/mskcc)
- [MSK-ACCESS](https://github.com/msk-access)
- [OncoKB](https://github.com/oncokb)
- [cBioPortal](https://github.com/cBioPortal)

## Contributors
 - See full list [here](https://github.com/ndbrown6/modules/graphs/contributors)
