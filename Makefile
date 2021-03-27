ifneq ("$(wildcard config.inc)", "")
	include config.inc
endif
ifneq ("$(wildcard project_config.inc)", "")
	include project_config.inc
endif
include innovation-lab/config/config.inc

export

NUM_ATTEMPTS ?= 10
NOW := $(shell date +"%F")
MAKELOG = log/$(@).$(NOW).log

USE_CLUSTER ?= true
QMAKE = innovation-lab/dodo-cloning-kit/runtime/qmake.pl -n $@.$(NOW) $(if $(SLACK_CHANNEL),-c $(SLACK_CHANNEL)) -r $(NUM_ATTEMPTS) -m -s -- make
NUM_JOBS ?= 100

define RUN_QMAKE
$(QMAKE) -e -f $1 -j $2 $(TARGET) && \
	mkdir -p completed_tasks && \
	touch completed_tasks/$@
endef

RUN_MAKE = $(if $(findstring false,$(USE_CLUSTER))$(findstring n,$(MAKEFLAGS)),+$(MAKE) -f $1,$(call RUN_QMAKE,$1,$(NUM_JOBS)))

#==================================================
# BETA testing
#==================================================

TARGETS += em_seq
em_seq :
	$(call RUN_MAKE,innovation-lab/workflows/em_seq.mk)
	
TARGETS += beadlink_prism
beadlink_prism :
	$(call RUN_MAKE,innovation-lab/workflows/beadlink_prism.mk)
	
TARGETS += fgbio_access
fgbio_access :
	$(call RUN_MAKE,innovation-lab/workflows/fgbio_access.mk)
	
TARGETS += pileup_metrics
pileup_metrics :
	$(call RUN_MAKE,innovation-lab/workflows/pileup_metrics.mk)
	
#==================================================
# MSK-ACCESS workflows
#==================================================

TARGETS += marianas_access
marianas_access :
	$(call RUN_MAKE,innovation-lab/workflows/marianas_access.mk)
	
TARGETS += genotype_access
genotype_access :
	$(call RUN_MAKE,innovation-lab/workflows/genotype_access.mk)
	
TARGETS += cluster_access
cluster_access :
	$(call RUN_MAKE,innovation-lab/workflows/cluster_access.mk)

#==================================================
# FASTQ / BAM file aligners
#==================================================

TARGETS += bwa_mem
bwa_mem :
	$(call RUN_MAKE,innovation-lab/aligners/bwa_mem.mk)

TARGETS += bismark_bt2
bismark_bt2 :
	$(call RUN_MAKE,innovation-lab/aligners/bismark_bt2.mk)

TARGETS += star_align
star_align :
	$(call RUN_MAKE,innovation-lab/aligners/star_align.mk)

#==================================================
# BAM file utilities
#==================================================

TARGETS += merge_alignments
merge_alignments :
	$(call RUN_MAKE,innovation-lab/bam_tools/merge_alignments.mk)

TARGETS += unmap_fasta
unmap_fasta :
	$(call RUN_MAKE,innovation-lab/bam_tools/unmap_fasta.mk)
	
TARGETS += unmap_bam
unmap_bam :
	$(call RUN_MAKE,innovation-lab/bam_tools/unmap_bam.mk)

TARGETS += merge_bam
merge_bam :
	$(call RUN_MAKE,innovation-lab/bam_tools/merge_bam.mk)
	
#==================================================
# FASTQ file utilities
#==================================================

TARGETS += extract_fastq
extract_fastq :
	$(call RUN_MAKE,innovation-lab/fastq_tools/extract_fastq.mk)

TARGETS += merge_fastq
merge_fastq :
	$(call RUN_MAKE,innovation-lab/fastq_tools/merge_fastq.mk)

TARGETS += subsample_fastq
subsample_fastq :
	$(call RUN_MAKE,innovation-lab/fastq_tools/subsample_fastq.mk)

#==================================================
# VCF file utilities
#==================================================

TARGETS += annotate_vcf_context
annotate_vcf_context :
	$(call RUN_MAKE,innovation-lab/vcf_tools/annotate_vcf_context.mk)
	
TARGETS += annotate_vcf_maf
annotate_vcf_maf :
	$(call RUN_MAKE,innovation-lab/vcf_tools/annotate_vcf_maf.mk)

#==================================================
# Copy number aberration callers
#==================================================

TARGETS += facets
facets :
	$(call RUN_MAKE,innovation-lab/copy_number/facets.mk)
	
TARGETS += ascat
ascat :
	$(call RUN_MAKE,innovation-lab/copy_number/ascat.mk)
	
TARGETS += cnvkit
cnvkit :
	$(call RUN_MAKE,innovation-lab/copy_number/cnvkit.mk)

#==================================================
# RNA structural variant/fusion callers
#==================================================

TARGETS += arriba
arriba :
	$(call RUN_MAKE,innovation-lab/structural_variants/arriba.mk)

TARGETS += defuse
defuse :
	$(call RUN_MAKE,innovation-lab/structural_variants/defuse.mk)
	
TARGETS += fusion_catcher
fusion_catcher :
	$(call RUN_MAKE,innovation-lab/structural_variants/fusioncatcher.mk)
	
TARGETS += star_fusion
star_fusion :
	$(call RUN_MAKE,innovation-lab/structural_variants/starfusion.mk)

#==================================================
# Quality control
#==================================================

TARGETS += rnaseq_metrics
rnaseq_metrics :
	$(call RUN_MAKE,innovation-lab/qc/rnaseq_metrics.mk)

TARGETS += fast_qc
fast_qc :
	$(call RUN_MAKE,innovation-lab/qc/fast_qc.mk)

TARGETS += cluster_samples
cluster_samples :
	$(call RUN_MAKE,innovation-lab/qc/cluster_samples.mk)

#==================================================
# RNA sequencing
#==================================================

TARGETS += umi_tools
umi_tools :
	$(call RUN_MAKE,innovation-lab/rna_seq/umi_tools.mk)
	
TARGETS += kallisto
kallisto :
	$(call RUN_MAKE,innovation-lab/rna_seq/kallisto.mk)

TARGETS += sum_reads
sum_reads :
	$(call RUN_MAKE,innovation-lab/rna_seq/sum_reads.mk)
	
TARGETS += sub_reads
sub_reads :
	$(call RUN_MAKE,innovation-lab/rna_seq/sub_reads.mk)

.PHONY : $(TARGETS)
