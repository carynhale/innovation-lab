ifndef CONFIG_INC

REF ?= b37
SEQ_PLATFORM ?= illumina
SAMPLE_SET_FILE ?= sample_sets.txt
SAMPLE_FILE ?= samples.txt
SAMPLE_SPLIT_FILE ?= samples.split.txt

# General interpreters
ANACONDA_PYTHON ?= python
PYTHON ?= python
PERL ?= perl
R ?= R
MY_RSCRIPT ?= Rscript
RSCRIPT ?= Rscript

# Python environment
ANACONDA_27_ENV ?= $(HOME)/share/usr/anaconda-envs/anaconda-2.7

# Anaconda environments
JRFLAB_MODULES_ENV ?= $(HOME)/share/usr/anaconda-envs/jrflab-modules-0.1.4
SUFAM_ENV ?= $(HOME)/share/usr/anaconda-envs/sufam-dev
MUTSIG_REPORT_ENV ?= $(HOME)/share/usr/anaconda-envs/mutsig-report-0.0.1
VEP_ENV ?= $(HOME)/share/usr/anaconda-envs/variant-effect-predictor-86
DECONSTRUCTSIGS_ENV ?= $(HOME)/share/usr/anaconda-envs/deconstructsigs
PHANGORN_ENV ?= $(HOME)/share/usr/anaconda-envs/phangorn
FGBIO_ENV ?= $(HOME)/share/usr/anaconda-envs/fgbio-0.8.1
HATCHET_ENV ?= $(HOME)/share/usr/anaconda-envs/hatchet
MACHINA_ENV ?= $(HOME)/share/usr/anaconda-envs/machina
CNTMD_ENV ?= $(HOME)/share/usr/anaconda-envs/cnt-md
CNTILP_ENV ?= $(HOME)/share/usr/anaconda-envs/cnt-ilp
PHYLO_ENV ?= $(HOME)/share/usr/anaconda-envs/phylotools
MEDICC_ENV ?= $(HOME)/share/usr/anaconda-envs/medicc
ASCAT_ENV ?= $(HOME)/share/usr/anaconda-envs/ascat

# General paths
JARDIR ?= $(HOME)/share/usr/lib/java
UNZIP ?= /usr/bin/unzip
JAVA_BIN ?=  $(JAVA8_BIN)
JAVA6_BIN ?= $(HOME)/share/usr/jdk1.6.0_45/bin/java
JAVA7_BIN ?= $(HOME)/share/usr/jdk1.7.0_45/bin/java
JAVA8_BIN ?= $(HOME)/share/usr/jdk1.8.0_121/bin/java

# Applications specific paths
FASTQC ?= $(PERL) $(HOME)/share/usr/FastQC/fastqc

MUTECT_JAR ?= $(JARDIR)/muTect-1.1.7.jar

BOWTIE ?= $(PERL) $(HOME)/share/usr/bowtie2-2.2.4/bowtie2
BWA ?= bwa

SAMTOOLS ?= samtools
SAMTOOLS2 ?= samtools

VCFUTILS ?= $(HOME)/share/usr/bin/vcfutils.pl

BCFTOOLS2 ?= bcftools
BCFTOOLS ?= bcftools

BEDTOOLS ?= $(HOME)/share/usr/bin/bedtools
BGZIP ?= $(HOME)/share/usr/bin/bgzip

IGVTOOLS ?= $(HOME)/share/usr/IGVTools/igvtools

VCFTOOLS ?= $(HOME)/share/usr/bin/vcftools-0.1.10

VCF_SORT ?= $(PERL) $(HOME)/share/usr/bin/vcfsorter.pl

SNP_EFF_JAR ?= $(JARDIR)/snpEff-4.3.jar
SNP_SIFT_JAR ?= $(JARDIR)/SnpSift-4.3.jar
SNP_EFF_CONFIG ?= modules/config/snpEff.conf

CUFFLINKS ?= cufflinks
CUFFCMP ?= cuffcompare
TOPHAT ?= tophat
DEFUSE ?= $(PERL) $(HOME)/share/usr/defuse-0.6.1/scripts/defuse.pl

ONCOFUSE_JAR ?= $(HOME)/share/usr/oncofuse-1.0.9b2/Oncofuse.jar
VARSCAN_JAR ?= $(JARDIR)/VarScan.v2.3.9.jar

PICARD_DIR ?= $(JARDIR)/picard-1.92
PICARD_JAR ?= $(JARDIR)/picard-tools-1.141/picard.jar

GATK_JAR ?= $(JARDIR)/GenomeAnalysisTK.jar
GATK_JAR2 ?= $(JARDIR)/GenomeAnalysisTK-3.7.jar

MARIANAS = /home/${USER}/share/usr/marianas-1.8.1/Marianas-1.8.1.jar
WALTZ = /home/${USER}/share/usr/marianas-1.8.1/Waltz-3.0.jar


# Resources
DB_NSFP ?= $(HOME)/share/reference/snpEff-4.1/dbNSFP3.0b1a.hg19.txt.gz

NSFP_FIELDS ?= Uniprot_acc_Polyphen2 Uniprot_id_Polyphen2 Polyphen2_HVAR_score Polyphen2_HVAR_pred 1000Gp3_AF ESP6500_AA_AF ESP6500_EA_AF MutationAssessor_pred MutationAssessor_score MutationTaster_pred MutationTaster_score PROVEAN_pred ExAC_Adj_AF clinvar_rs clinvar_clnsig Interpro_domain




# scripts
SCRIPTS_DIR ?= modules/scripts
MERGE ?= $(SCRIPTS_DIR)/merge.R
VCF_TO_TABLE ?= $(SCRIPTS_DIR)/vcfToTable.R
INTRON_POSN_LOOKUP ?= $(SCRIPTS_DIR)/posnGeneLookup.pl
RBIND ?= $(SCRIPTS_DIR)/rbind.R
NORMAL_FILTER ?= $(PERL) $(SCRIPTS_DIR)/normalFilterVCF.pl
SOMATIC_FILTER_VCF ?= $(PERL) $(SCRIPTS_DIR)/somaticFilterVCF.pl

GET_INSERT_SIZE ?= $(HOME)/share/usr/bin/getInsertSize.py


VT ?= $(HOME)/share/usr/bin/vt
TABIX ?= $(HOME)/share/usr/bin/tabix

CLUSTER_ENGINE ?= SGE
ifeq ($(CLUSTER_ENGINE),PBS)
MYSQLDB_HOST = cpu-6-2
EMBL_MYSQLDB_HOST = cpu-6-2
EMBL_MYSQLDB_PORT = 9992
EMBL_MYSQLDB_USER = embl
EMBL_MYSQLDB_PW = embl
EMBL_MYSQLDB_DB = homo_sapiens_core_85_37
endif
ifeq ($(CLUSTER_ENGINE),LSF)
MYSQLDB_HOST = swan
EMBL_MYSQLDB_HOST = swan
EMBL_MYSQLDB_PORT = 38493
EMBL_MYSQLDB_USER = embl
EMBL_MYSQLDB_PW = embl
EMBL_MYSQLDB_DB = homo_sapiens_core_75_37
endif
ifeq ($(CLUSTER_ENGINE),SGE)
MYSQLDB_HOST = ika
EMBL_MYSQLDB_HOST = ika
EMBL_MYSQLDB_PORT = 38493
EMBL_MYSQLDB_USER = embl
EMBL_MYSQLDB_PW = embl
EMBL_MYSQLDB_DB = homo_sapiens_core_75_37
endif
SNVBOX_CONF = $(PWD)/modules/db/snv_box-$(MYSQLDB_HOST).conf
FATHMM_CONF = $(PWD)/modules/db/fathmm_config-$(MYSQLDB_HOST).ini

INTEGRATE ?= $(HOME)/share/usr/Integrate_0.2.0/bin/Integrate

VCF2MAF ?= $(PERL) $(HOME)/share/usr/bin/vcf2maf.pl
MAF2MAF ?= $(PERL) $(HOME)/share/usr/bin/maf2maf.pl
MAF2VCF ?= $(PERL) $(HOME)/share/usr/bin/maf2vcf.pl

HISAT ?= $(HOME)/share/usr/hisat2-2.0.2-beta/hisat2

TVC ?= $(HOME)/share/usr/bin/tvc

ANNOVAR = $(PERL) $(HOME)/share/usr/annovar-2017-07-16/table_annovar.pl

VEP_PATH = $(VEP_ENV)/bin

SPLIT_BED = python modules/scripts/split_bed.py

INDEL_FILTER_VCF = python modules/vcf_tools/indel_filter_vcf.py
SNP_FILTER_VCF = python modules/vcf_tools/snp_filter_vcf.py

MERGE_VCF = python modules/vcf_tools/merge_vcf.py

MEDICC_VAR = modules/config/medicc.conf
MEDICC_BIN = $(MEDICC_ENV)/opt/medicc

CNTILP_CONF = modules/config/cnt-ilp.conf
CNTILP_BIN = $(CNTILP_ENV)/opt/CNT-ILP/build

CNTMD_CONF = modules/config/cnt-md.conf
CNTMD_BIN = $(CNTMD_ENV)/opt/CNT-MD/build

MACHINA_VAR = modules/config/machina.conf
MACHINA_BIN = $(MACHINA_ENV)/opt/machina/build

HATCHET_VAR = modules/config/hatchet.conf
HATCHET_BIN = $(HATCHET_ENV)/opt/machina/build

MARIANAS_UMI_LENGTH = 3
MARIANAS_MIN_MAPQ = 1
MARIANAS_MIN_BAQ = 20
MARIANAS_MISMATCH = 0
MARIANAS_WOBBLE = 1
MARIANAS_MIN_CONSENSUS = 90

WALTZ_MIN_MAPQ = 20
WALTZ_BED_FILE = $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.waltz.bed

endif
CONFIG_INC = true
