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
JAVA ?=  $(JAVA8)

# java applications
MUTECT_JAR ?= $(JARDIR)/muTect-1.1.7.jar
VARSCAN_JAR ?= $(JARDIR)/VarScan.v2.3.9.jar
MARIANAS ?= $(HOME)/share/usr/marianas-1.8.1/Marianas-1.8.1.jar
WALTZ ?= $(HOME)/share/usr/marianas-1.8.1/Waltz-3.0.jar
PICARD_DIR ?= $(JARDIR)/picard-1.92

# SCRIPTS_DIR ?= modules/dodo-cloning-kit
# VCF2MAF ?= $(PERL) $(HOME)/share/usr/bin/vcf2maf.pl
# MAF2MAF ?= $(PERL) $(HOME)/share/usr/bin/maf2maf.pl
# MAF2VCF ?= $(PERL) $(HOME)/share/usr/bin/maf2vcf.pl

REF ?= b37
SEQ_PLATFORM ?= illumina

endif
CONFIG_INC = true
