# innovation-lab
[![Build Status]](https://travis-ci.com/ndbrown6/innovation-lab/)

## Introduction
This is a lightweight fork of the jrflab [modules](https://github.com/jrflab/modules).

## Installation
The easiest way to download this pipeline is to clone the repository.

```
git clone --recurse-submodules https://github.com/ndbrown6/modules.git
```

## Dependencies
- An instance of [anaconda](https://www.anaconda.com) or [miniconda](https://conda.io/en/latest/miniconda.html)
- IMB's Platform Load Sharing Facility (LSF) for resource management
- [Full dependencies](https://github.com/ndbrown6/modules/tree/master/conda)

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
- Mutation calling:
	* [MuTect](https://www.nature.com/articles/nbt.2514)
	* [Platypus](https://www.nature.com/articles/ng.3036)
	* [Strelka](https://academic.oup.com/bioinformatics/article/28/14/1811/218573)
	* [Scalpel](https://doi.org/10.1038/nprot.2016.150)
	* [Lancet](https://www.nature.com/articles/s42003-018-0023-9)
	* [Varscan](https://academic.oup.com/bioinformatics/article/25/17/2283/210190)
- Copy number aberrations:
	* [FACETS](https://academic.oup.com/nar/article/44/16/e131/2460163)
	* [ASCAT](https://www.pnas.org/content/107/39/16910)
	* [CNVkit](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004873)
- Annotation:
	* Annovar
	* snpEff
	* SIFT
	* pph2
	* [vcf2maf](https://github.com/mskcc/vcf2maf)
	* VEP
	* [OncoKB](https://github.com/oncokb/oncokb-annotator)
	* ClinVar
- Miscellaneous:
	* [HLA Polysolver](https://www.nature.com/articles/nbt.3344)
	* [MSIsensor](https://academic.oup.com/bioinformatics/article/30/7/1015/236553)

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
- [MSK-ACCESS](https://github.com/msk-access)
- [OncoKB](https://github.com/oncokb)
- [cBioPortal](https://github.com/cBioPortal)
- [GRAIL](https://github.com/grailbio)

## Contributors
 - See full list [here](https://github.com/ndbrown6/modules/graphs/contributors)
 
## Honorable mentions
- With contributions from Nadeem Riaz & Fong Chun Chan

