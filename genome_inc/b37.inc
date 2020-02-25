ifndef B37_INC
REF = b37

REF_DIR = $(HOME)/share/lib/ref_files/b37
REF_FASTA ?= $(REF_DIR)/human_g1k_v37.fasta
DMP_FASTA ?= $(REF_DIR)/human_g1k_v37_dmp.fasta
DECOY_FASTA ?= $(REF_DIR)/human_g1k_v37_decoy.fasta
DBSNP ?= $(REF_DIR)/dbsnp_138.b37.vcf
KNOWN_INDELS ?= $(REF_DIR)/Mills_and_1000G_gold_standard.indels.b37.vcf
OMNI ?= $(REF_DIR)/1000G_omni2.5.b37.vcf
HAPMAP ?= $(REF_DIR)/hapmap_3.3.b37.vcf
HAPMAP_POP_FILE ?= $(REF_DIR)/hapmap_3.3_b37_pop_stratified_af.vcf.gz
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE) | cut -f1 | sort | uniq)
REF_MSI = $(REF_DIR)/human_g1k_v37_msi.list
ENSEMBLE = $(HOME)/share/lib/resource_files/Ensemble_b37.RData

endif
B37_INC = true
