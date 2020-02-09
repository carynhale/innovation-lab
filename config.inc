ifndef CONFIG_INC

# conda
INNOVATION_LAB_001 ?= $(HOME)/share/usr/env/innovation-lab-0.0.1

# alias
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
PICRAD ?= picard
GATK ?= gatk3

# java
JAVA6 ?= $(HOME)/share/usr/java/jdk1.6.0_45/bin/java
JAVA7 ?= $(HOME)/share/usr/java/jdk1.7.0_80/bin/java
JAVA8 ?= $(HOME)/share/usr/java/jdk1.8.0_202/bin/java

# uncategorized
REF ?= b37
SEQ_PLATFORM ?= illumina
SCRIPTS ?= modules/dodo-cloning-kit

endif
CONFIG_INC = true
