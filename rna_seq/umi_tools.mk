include innovation-lab/Makefile.inc
include innovation-lab/config/align.inc
include innovation-lab/config/gatk.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR = log/umi_tools.$(NOW)

ALIGNER := star
BAM_NO_REALN = true
BAM_NO_RECAL = true
BAM_NO_FILTER = true
BAM_DUP_TYPE = none
SEQ_PLATFROM = illumina

UMI_PATTERN = NNNX

STAR_OPTS = --genomeDir $(STAR_REF) \
            --outSAMtype BAM SortedByCoordinate \
			--twopassMode Basic \
            --outReadsUnmapped None \
            --outSAMunmapped Within \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 12 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 200000 \
            --alignIntronMax 200000 \
            --chimSegmentReadGapMax parameter 3 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --chimOutType WithinBAM \
			--quantMode GeneCounts

umi_tools : $(foreach sample,$(SAMPLES),umi_tools/$(sample)/$(sample)_R1.fastq.gz) \
			$(foreach sample,$(SAMPLES),umi_tools/$(sample)/$(sample)_R1_cl.fastq.gz) \
			$(foreach sample,$(SAMPLES),umi_tools/$(sample)/$(sample)_R2_cl.fastq.gz) \
			$(foreach sample,$(SAMPLES),star/$(sample).Aligned.sortedByCoord.out.bam) \
			$(foreach sample,$(SAMPLES),star/$(sample).Aligned.sortedByCoord.out.bam.bai) \
			$(foreach sample,$(SAMPLES),bam/$(sample).bam)

define copy-fastq
umi_tools/$1/$1_R1.fastq.gz : $3
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
									 mkdir -p umi_tools/$1 && \
								     $(RSCRIPT) $(SCRIPTS_DIR)/fastq_tools/copy_fastq.R \
								     --sample_name $1 \
								     --directory_name umi_tools \
								     --fastq_files '$$^'")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call copy-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))
	
define clip-fastq
umi_tools/$1/$1_R1_cl.fastq.gz : umi_tools/$1/$1_R1.fastq.gz
	$$(call RUN,-c -n 1 -s 8G -m 16G -v $(UMITOOLS_ENV),"set -o pipefail && \
														 umi_tools extract \
														 --bc-pattern=$$(UMI_PATTERN) \
														 -I $$(<) \
														 --stdout $$(@) \
														 --log=umi_tools/$1/$1_R1_cl.log")

umi_tools/$1/$1_R2_cl.fastq.gz : umi_tools/$1/$1_R2.fastq.gz
	$$(call RUN,-c -n 1 -s 8G -m 16G -v $(UMITOOLS_ENV),"set -o pipefail && \
														 umi_tools extract \
														 --bc-pattern=$$(UMI_PATTERN) \
														 -I $$(<) \
														 --stdout $$(@) \
														 --log=umi_tools/$1/$1_R2_cl.log")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call clip-fastq,$(sample))))

define align-fastq
star/$1.Aligned.sortedByCoord.out.bam : umi_tools/$1/$1_R1_cl.fastq.gz umi_tools/$1/$1_R2_cl.fastq.gz
	$$(call RUN,-n 4 -s 6G -m 10G -w 1440,"set -o pipefail && \
										   STAR $$(STAR_OPTS) \
										   --outFileNamePrefix star/$1. \
										   --runThreadN 4 \
										   --outSAMattrRGline \"ID:$1\" \"LB:$1\" \"SM:$1\" \"PL:$${SEQ_PLATFORM}\" \
										   --readFilesIn $$(<) $$(<<) \
										   --readFilesCommand zcat")
                                   
star/$1.Aligned.sortedByCoord.out.bam.bai : star/$1.Aligned.sortedByCoord.out.bam
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
                                  $$(SAMTOOLS) index $$(<)")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call align-fastq,$(sample))))

define dedup-bam
bam/$1.bam : star/$1.Aligned.sortedByCoord.out.bam
	$$(call RUN,-c -n 1 -s 24G -m 48G -v $(UMITOOLS_ENV) -w 1440,"set -o pipefail && \
																  umi_tools dedup \
																  -I $$(<) \
																  -S $$(@) \
																  --output-stats=umi_tools/$1/$1 \
																  --paired")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call dedup-bam,$(sample))))
    
..DUMMY := $(shell mkdir -p version; \
			 $(UMITOOLS_ENV)/bin/umi_tools --version > version/umi_tools.txt; \
             echo "STAR" >> version/umi_tools.txt; \
			 STAR --version >> version/umi_tools.txt; \
             $(SAMTOOLS) --version >> version/umi_tools.txt)
.SECONDARY: 
.DELETE_ON_ERROR:
.PHONY: umi_tools

include innovation-lab/fastq_tools/fastq.mk
include innovation-lab/bam_tools/process_bam.mk
include innovation-lab/aligners/align.mk
