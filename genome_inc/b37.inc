ifndef B37_INC
REF = b37

REF_DIR = $(HOME)/share/lib/ref_files
REF_FASTA ?= $(REF_DIR)/.../human_g1k_v37.fasta
DMP_FASTA ?= $(REF_DIR)/.../b37.fasta
NOCHR_REF_FASTA ?= $(REF_DIR)/.../human_g1k_v37.fasta
CHR_REF_FASTA ?= $(REF_DIR)/.../ucsc.hg19.fasta
REF_2BIT ?= $(REF_DIR)/.../human_g1k_v37.2bit
DBSNP ?= $(REF_DIR)/.../dbsnp_138.b37.gmaf.vcf.gz
KNOWN_INDELS ?= $(REF_DIR)/.../Mills_and_1000G_gold_standard.indels.b37.vcf.gz
OMNI ?= $(REF_DIR)/.../1000G_omni2.5.b37.vcf.gz
HAPMAP ?= $(REF_DIR)/.../hapmap_3.3.b37.vcf.gz
COSMIC ?= $(REF_DIR)/.../cosmic_v82/CosmicCodingMuts.vcf.gz
COSMIC_NONCODING ?= $(REF_DIR)/.../cosmic_v82/CosmicNonCodingVariants.vcf.gz
HAPMAP_POP_FILE ?= $(REF_DIR)/.../hapmap_3.3_b37_pop_stratified_af.vcf.gz

CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE) | cut -f1 | sort | uniq)

REF_MSI = $(REF_DIR)/.../b37_msi.list

VEP_REF = GRCh37
VEP_DATA = $(REF_DIR)/vep/v86
VEP_OVERRIDES = $(REF_DIR)/isoform_overrides_at_mskcc
VEP_SPECIES = homo_sapiens
ONCOKB = $(REF_DIR)/allActionableVariants.txt

endif
B37_INC = true

