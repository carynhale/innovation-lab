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

# java command lines
JAVA6 ?= $(HOME)/share/usr/java/jdk1.6.0_45/bin/java -Djava.io.tmpdir=$(TMPDIR)
JAVA7 ?= $(HOME)/share/usr/java/jdk1.7.0_80/bin/java -Djava.io.tmpdir=$(TMPDIR)
JAVA8 ?= $(HOME)/share/usr/java/jdk1.8.0_202/bin/java -Djava.io.tmpdir=$(TMPDIR)

# picard command lines
SORT_SAM = $(JAVA) -Xmx$1 -jar $(PICARD_DIR)/SortSam.jar $(PICARD_OPTS) TMP_DIR=$(TMPDIR)
REORDER_SAM = $(JAVA) -Xmx$(1) -jar $(PICARD_DIR)/ReorderSam.jar $(PICARD_OPTS)
MARK_DUP = $(JAVA) -Xmx$(1) -jar $(PICARD_DIR)/MarkDuplicates.jar $(PICARD_OPTS)
MERGE_SAM = $(JAVA) -Xmx$(PICARD_MEM) -jar $(PICARD_DIR)/MergeSamFiles.jar $(PICARD_OPTS)
INTERVAL_LIST_TOOL = $(JAVA) -Xmx$(PICARD_MEM) -jar $(PICARD_DIR)/IntervalListTools.jar $(PICARD_OPTS)
CREATE_SEQ_DICT = $(JAVA) -Xmx$(PICARD_MEM) -jar $(PICARD_DIR)/CreateSequenceDictionary.jar $(PICARD_OPTS)
CALC_HS_METRICS = $(JAVA) -Xmx$(PICARD_MEM) -jar $(PICARD_DIR)/CalculateHsMetrics.jar $(PICARD_OPTS)
COLLECT_MULT_METRICS = $(JAVA) -Xmx$(PICARD_MEM) -jar $(PICARD_DIR)/CollectMultipleMetrics.jar $(PICARD_OPTS)
COLLECT_TARGETED_METRICS = $(JAVA) -Xmx$(PICARD_MEM) -jar $(PICARD_DIR)/CollectTargetedPcrMetrics.jar $(PICARD_OPTS)
FIX_MATE = $(JAVA) -Xmx$(1) -jar $(PICARD_DIR)/FixMateInformation.jar $(PICARD_OPTS) TMP_DIR=$(TMPDIR)
SAM_TO_FASTQ = $(JAVA) -Xmx$(1) -jar $(PICARD_DIR)/SamToFastq.jar $(PICARD_OPTS)
CLEAN_BAM = $(JAVA) -Xmx$(1) -jar $(PICARD_DIR)/CleanSam.jar $(PICARD_OPTS)
ADD_RG_MEM = $(JAVA) -Xmx$(1) -jar $(PICARD_DIR)/AddOrReplaceReadGroups.jar $(PICARD_OPTS)
VALIDATE_SAM = $(JAVA) -Xmx$(PICARD_MEM) -jar $(PICARD_DIR)/ValidateSamFile.jar $(PICARD_OPTS)

# gatk command lines
ANALYZE_COV = GATK -Xmx10G -T AnalyzeCovariates 

# paths
SCRIPTS_DIR ?= innovation-lab/dodo-cloning-kit

endif
CONFIG_INC = true
