include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/msk_access.$(NOW)

# MSK_ACCESS_WORKFLOW += align_fastq
# MSK_ACCESS_WORKFLOW += umi_collapse
# MSK_ACCESS_WORKFLOW += align_collapsed
# MSK_ACCESS_WORKFLOW += copy_bam
# MSK_ACCESS_WORKFLOW += interval_metrics
# MSK_ACCESS_WORKFLOW += umi_qc
# MSK_ACCESS_WORKFLOW += plot_metrics
# MSK_ACCESS_WORKFLOW += cluster_samples

msk_access : $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)_R1.fastq.gz) \
		   	 $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)_R2.fastq.gz) \
		   	 $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)_R1_umi-clipped.fastq.gz) \
		   	 $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)_R2_umi-clipped.fastq.gz)

MARIANAS_UMI_LENGTH ?= 3
MARIANAS_MIN_MAPQ ?= 1
MARIANAS_MIN_BAQ ?= 20
MARIANAS_MISMATCH ?= 0
MARIANAS_WOBBLE ?= 1
MARIANAS_MIN_CONSENSUS ?= 90
WALTZ_MIN_MAPQ ?= 20

define copy-fastq
marianas/$1/$1_R1.fastq.gz marianas/$1/$1_R2.fastq.gz : $3
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
									 mkdir -p marianas/$1 && \
								     $(RSCRIPT) $(SCRIPTS_DIR)/fastq_tools/copyfastq.R --sample_name $1 --fastq_files '$$^'")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call copy-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))
 

define clip-umi
marianas/$1/$1_R1_umi-clipped.fastq.gz marianas/$1/$1_R2_umi-clipped.fastq.gz : marianas/$1/$1_R1.fastq.gz marianas/$1/$1_R2.fastq.gz
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  cd marianas/$1/ && \
									  $(call MARIANAS_CMD,2G,8G) \
									  org.mskcc.marianas.umi.duplex.fastqprocessing.ProcessLoopUMIFastq \
									  $1_R1.fastq.gz $1_R2.fastq.gz \
									  $(MARIANAS_UMI_LENGTH) && \
									  cd ../..")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call clip-umi,$(sample))))

# include modules/test/bam_tools/alignfastq.mk
# include modules/test/bam_tools/collapseumi.mk
# include modules/test/bam_tools/aligncollapsed.mk
# include modules/test/bam_tools/copybam.mk
# include modules/test/qc/intervalmetrics.mk
# include modules/test/qc/umiqc.mk
# include modules/test/qc/plotmetrics.mk
# include modules/test/qc/clustersamples.mk

..DUMMY := $(shell mkdir -p version)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: msk_access
