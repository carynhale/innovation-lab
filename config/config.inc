ifndef CONFIG_INC

REF ?= b37
SEQ_PLATFORM ?= illumina

# conda
INNOVATION_LAB_ENV ?= $(HOME)/share/usr/env/innovation-lab-0.0.1
SUPERHEAT_ENV = $(HOME)/share/usr/env/r-complexheatmap-2.2.0
VARIANT_ANNOTATION_ENV = $(HOME)/share/usr/env/r-variantannotation-1.32.0
ANNOTATION_ENV = $(HOME)/share/usr/env/r-bsgenome.hsapiens.ucsc.hg19-1.4.0
SUMREADS_ENV = $(HOME)/share/usr/env/r-bsgenome.hsapiens.ucsc.hg19-1.4.0
GENOTYPE_VARIANTS_ENV ?= $(HOME)/share/usr/env/genotype-variants
VCF2MAF_ENV ?= $(HOME)/share/usr/env/vcf2maf-1.6.17
DEFUSE_ENV ?= $(HOME)/share/usr/env/defuse-0.8.0
FUSIONCATCHER_ENV ?= $(HOME)/share/usr/env/fusioncatcher-1.2.0
STARFUSION_ENV ?= $(HOME)/share/usr/env/starfusion-1.6.0
ABSOLUTE_ENV ?= $(HOME)/share/usr/env/cntu-0.0.1
BISMARK_ENV ?= $(HOME)/share/usr/env/bismark-0.22.1
UMITOOLS_ENV ?= $(HOME)/share/usr/env/umi-tools-1.0.1
KALLISTO_ENV ?= $(HOME)/share/usr/env/kallisto-0.46.2
ARRIBA_ENV ?= $(HOME)/share/usr/env/arriba-2.0.0
FASTQC_ENV ?= $(HOME)/share/usr/env/cutadapt-3.2
SEQTK_ENV ?= $(HOME)/share/usr/env/seqtk-1.3
SALMON_ENV ?= $(HOME)/share/usr/env/salmon-1.4.0
CNVKIT_ENV ?= $(HOME)/share/usr/env/cnvkit-0.9.8
MSISENSOR_ENV ?= $(HOME)/share/usr/env/msisensor-pro-1.1.a
GATK_ENV ?= $(HOME)/share/usr/env/gatk-3.8-10
SUFAM_ENV ?= $(HOME)/share/usr/env/sufam-v0.4.2
FASTQ_SPLITTER_ENV ?= $(HOME)/share/usr/env/fastqsplitter-1.2.0
FACETS_ENV ?= $(HOME)/share/usr/env/r-facets-0.6.2
SAMBAMBA_ENV ?= $(HOME)/share/usr/env/sambamba-0.8.2
CRAVAT_ENV ?= $(HOME)/share/usr/env/open-cravat-2.2.7

# aliases
PYTHON ?= python
PERL ?= perl
RSCRIPT ?= Rscript
R ?= R
BWA ?= bwa
SAMTOOLS ?= samtools
SAMBAMBA ?= sambamba
BEDTOOLS ?= bedtools
BGZIP ?= bgzip
UNZIP ?= unzip
TABIX ?= tabix
PICARD ?= picard
GATK ?= gatk
BAM2FASTX ?= bam2fastx
GENOTYPE_VARIANTS ?= genotype_variants
VCF2MAF ?= vcf2maf.pl
MAF2VCF ?= maf2vcf.pl
DEFUSE ?= defuse_run.pl
STAR_FUSION ?= STAR-Fusion
FASTQC ?= fastqc
CUTADAPT ?= cutadapt
SEQTK ?= seqtk
SUFAM ?= sufam


# java command lines
JAVA6 ?= $(HOME)/share/usr/java/jdk1.6.0_45/bin/java -Djava.io.tmpdir=$(TMPDIR)
JAVA7 ?= $(HOME)/share/usr/java/jdk1.7.0_80/bin/java -Djava.io.tmpdir=$(TMPDIR)
JAVA8 ?= $(HOME)/share/usr/java/jdk1.8.0_74/bin/java -Djava.io.tmpdir=$(TMPDIR)

# picard command lines
PICARD_MEM = 16G
PICARD_OPTS = VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000 TMP_DIR=$(TMPDIR)
SORT_SAM = $(PICARD) -Xmx$(PICARD_MEM) SortSam $(PICARD_OPTS)
REORDER_SAM = $(PICARD) -Xmx$(PICARD_MEM) ReorderSam $(PICARD_OPTS)
MARK_ADAPTERS = $(PICARD) -Xmx$(PICARD_MEM) MarkIlluminaAdapters $(PICARD_OPTS)
MARK_DUP = $(PICARD) -Xmx$(PICARD_MEM) MarkDuplicates $(PICARD_OPTS)
MERGE_SAM = $(PICARD) -Xmx$(PICARD_MEM) MergeSamFiles $(PICARD_OPTS)
INTERVAL_LIST_TOOL = $(PICARD) -Xmx$(PICARD_MEM) IntervalListTools $(PICARD_OPTS)
CREATE_SEQ_DICT = $(PICARD) -Xmx$(PICARD_MEM) CreateSequenceDictionary $(PICARD_OPTS)
CALC_HS_METRICS = $(PICARD) -Xmx$(PICARD_MEM) CollectHsMetrics $(PICARD_OPTS)
COLLECT_MULT_METRICS = $(PICARD) -Xmx$(PICARD_MEM) CollectMultipleMetrics $(PICARD_OPTS)
COLLECT_TARGETED_METRICS = $(PICARD) -Xmx$(PICARD_MEM) CollectTargetedPcrMetrics $(PICARD_OPTS)
COLLECT_RNASEQ_METRICS = $(PICARD) -Xmx$(PICARD_MEM) CollectRnaSeqMetrics $(PICARD_OPTS)
COLLECT_RRBS_METRICS = $(PICARD) -Xmx$(PICARD_MEM) CollectRrbsMetrics $(PICARD_OPTS)
COLLECT_ALIGNMENT_METRICS = $(PICARD) -Xmx$(PICARD_MEM) CollectAlignmentSummaryMetrics $(PICAD_OPTS)
COLLECT_INSERT_METRICS = $(PICARD) -Xmx$(PICARD_MEM) CollectInsertSizeMetrics $(PICAD_OPTS)
COLLECT_OXOG_METRICS = $(PICARD) -Xmx$(PICARD_MEM) CollectOxoGMetrics $(PICAD_OPTS)
COLLECT_GC_BIAS = $(PICARD) -Xmx$(PICARD_MEM) CollectGcBiasMetrics $(PICAD_OPTS)
CALC_LIBRARY_COMPLEXITY = $(PICARD) -Xmx$(PICARD_MEM) EstimateLibraryComplexity $(PICAD_OPTS)
FIX_MATE = $(PICARD) -Xmx$(PICARD_MEM) FixMateInformation $(PICARD_OPTS)
SAM_TO_FASTQ = $(PICARD) -Xmx$(PICARD_MEM) SamToFastq $(PICARD_OPTS)
FASTQ_TO_SAM = $(PICARD) -Xmx$(PICARD_MEM) FastqToSam $(PICARD_OPTS)
CLEAN_BAM = $(PICARD) -Xmx$(PICARD_MEM) CleanSam $(PICARD_OPTS)
ADD_RG = $(PICARD) -Xmx$(PICARD_MEM) AddOrReplaceReadGroups $(PICARD_OPTS)
VALIDATE_SAM = $(PICARD) -Xmx$(PICARD_MEM) ValidateSamFile $(PICARD_OPTS)
REVERT_SAM = $(PICARD) -Xmx$(PICARD_MEM) RevertSam $(PICARD_OPTS)
MERGE_ALIGNMENTS = $(PICARD) -Xmx$(PICARD_MEM) MergeBamAlignment $(PICAD_OPTS)
BAM_INDEX = $(PICARD) -Xmx$(PICARD_MEM) BamIndexStats $(PICAD_OPTS)

# gatk command lines
GATK_MEM = 16G
GATK_OPTS = -S LENIENT
ANALYZE_COV = $(GATK) -Xmx$(GATK_MEM) -T AnalyzeCovariates $(GATK_OPT)

# gbc command line
GBC ?= $(HOME)/share/usr/GetBaseCounts/GetBaseCounts

# waltz command line
WALTZ ?= $(HOME)/share/usr/marianas-1.8.1/Waltz-3.0.jar

# fgbio command line
FGBIO ?= $(HOME)/share/usr/fgbio-1.1.0/fgbio-1.1.0.jar

# mutect command line
MUTECT_EXE = $(HOME)/share/usr/mutect-1.1.7/muTect-1.1.7.jar
MUTECT ?= $(JAVA7) -Xmx11G -jar $(MUTECT_EXE) --analysis_type MuTect

# paths
SCRIPTS_DIR ?= innovation-lab/dodo-cloning-kit
TMPDIR ?= $(HOME)/share/tmp

endif
CONFIG_INC = true
