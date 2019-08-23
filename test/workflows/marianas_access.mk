include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/marianas_access.$(NOW)
PHONY += marianas bam metrics

MSK_ACCESS_WORKFLOW += clip_umi
MSK_ACCESS_WORKFLOW += align_fastq
MSK_ACCESS_WORKFLOW += umi_collapse
MSK_ACCESS_WORKFLOW += align_collapsed
MSK_ACCESS_WORKFLOW += copy_bam
MSK_ACCESS_WORKFLOW += interval_metrics
#MSK_ACCESS_WORKFLOW += umi_qc

msk_access_workflow : $(MSK_ACCESS_WORKFLOW)

include modules/test/fastq_tools/clipumi.mk
include modules/test/bam_tools/alignfastq.mk
include modules/test/bam_tools/collapseumi.mk
include modules/test/bam_tools/aligncollapsed.mk
include modules/test/bam_tools/copybam.mk
include modules/test/qc/intervalmetrics.mk
#include modules/test/qc/umiqc.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
