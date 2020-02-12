ifndef CONFIG_INC

REF ?= b37
SEQ_PLATFORM ?= illumina

# conda
INNOVATION_LAB_ENV ?= $(HOME)/share/usr/env/innovation-lab-0.0.1

# aliases
PYTHON ?= python
PERL ?= perl
RSCRIPT ?= Rscript
R ?= R
BWA ?= bwa
SAMTOOLS ?= samtools
BEDTOOLS ?= bedtools
BGZIP ?= bgzip
UNZIP ?= unzip
TABIX ?= tabix
PICARD ?= picard
GATK ?= gatk3
BAM2FASTX ?= bam2fastx
GETBASECOUNTS ?= $(HOME)/share/usr/bin/GetBaseCounts

# java command lines
JAVA6 ?= $(HOME)/share/usr/java/jdk1.6.0_45/bin/java -Djava.io.tmpdir=$(TMPDIR)
JAVA7 ?= $(HOME)/share/usr/java/jdk1.7.0_80/bin/java -Djava.io.tmpdir=$(TMPDIR)
JAVA8 ?= $(HOME)/share/usr/java/jdk1.8.0_74/bin/java -Djava.io.tmpdir=$(TMPDIR)

# picard command lines
PICARD_MEM = 12G
PICARD_OPTS = VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000 TMP_DIR=$(TMPDIR)
SORT_SAM = $(PICARD) -Xmx$(PICARD_MEM) SortSam $(PICARD_OPTS)
REORDER_SAM = $(PICARD) -Xmx$(PICARD_MEM) ReorderSam $(PICARD_OPTS)
MARK_DUP = $(PICARD) -Xmx$(PICARD_MEM) MarkDuplicates $(PICARD_OPTS)
MERGE_SAM = $(PICARD) -Xmx$(PICARD_MEM) MergeSamFiles $(PICARD_OPTS)
INTERVAL_LIST_TOOL = $(PICARD) -Xmx$(PICARD_MEM) IntervalListTools $(PICARD_OPTS)
CREATE_SEQ_DICT = $(PICARD) -Xmx$(PICARD_MEM) CreateSequenceDictionary $(PICARD_OPTS)
CALC_HS_METRICS = $(PICARD) -Xmx$(PICARD_MEM) CollectHsMetrics $(PICARD_OPTS)
COLLECT_MULT_METRICS = $(PICARD) -Xmx$(PICARD_MEM) CollectMultipleMetrics $(PICARD_OPTS)
COLLECT_TARGETED_METRICS = $(PICARD) -Xmx$(PICARD_MEM) CollectTargetedPcrMetrics $(PICARD_OPTS)
FIX_MATE = $(PICARD) -Xmx$(PICARD_MEM) FixMateInformation $(PICARD_OPTS)
SAM_TO_FASTQ = $(PICARD) -Xmx$(PICARD_MEM) SamToFastq $(PICARD_OPTS)
CLEAN_BAM = $(PICARD) -Xmx$(PICARD_MEM) CleanSam $(PICARD_OPTS)
ADD_RG = $(PICARD) -Xmx$(PICARD_MEM) AddOrReplaceReadGroups $(PICARD_OPTS)
VALIDATE_SAM = $(PICARD) -Xmx$(PICARD_MEM) ValidateSamFile $(PICARD_OPTS)
REVERT_SAM = $(PICARD) -Xmx$(PICARD_MEM) RevertSam $(PICARD_OPTS)
MERGE_ALIGNMENTS = $(PICARD) -Xmx$(PICARD_MEM) MergeBamAlignment $(PICAD_OPTS)

# gatk command lines
GATK_MEM = 18G
GATK_OPTS = -S LENIENT
ANALYZE_COV = $(GATK) -Xmx$(GATK_MEM) -T AnalyzeCovariates $(GATK_OPT)

# marianas command line
MARIANAS ?= $(HOME)/share/usr/marianas-1.8.1/Marianas-1.8.1.jar

# waltz command line
WALTZ ?= $(HOME)/share/usr/marianas-1.8.1/Waltz-3.0.jar

# paths
SCRIPTS_DIR ?= innovation-lab/dodo-cloning-kit
TMPDIR ?= $(HOME)/share/tmp

endif
CONFIG_INC = true
