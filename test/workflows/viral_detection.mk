include modules/Makefile.inc
include modules/config.inc

LOGDIR = log/viral_detection.$(NOW)
PHONY += unmapped_reads

VIRUS_WORKFLOW += extract_unmapped
VIRUS_WORKFLOW += bam_to_fasta
VIRUS_WORKFLOW += blast_reads
VIRUS_WORKFLOW += krona_classify

viral_detection_workflow : $(VIRUS_WORKFLOW)

include modules/bam_tools/extract_unmapped.mk
include modules/bam_tools/bam_to_fasta.mk
include modules/aligners/blast_aligner.mk
include modules/virus/krona_classify.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)


