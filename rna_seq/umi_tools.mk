include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR = log/umi_tools.$(NOW)

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
	    $(foreach sample,$(SAMPLES),umi_tools/$(sample)/$(sample)_R2.fastq.gz) \
	    $(foreach sample,$(SAMPLES),umi_tools/$(sample)/$(sample)_R1_cl.fastq.gz) \
	    $(foreach sample,$(SAMPLES),star/$(sample).Aligned.sortedByCoord.out.bam) \
	    $(foreach sample,$(SAMPLES),star/$(sample).Aligned.sortedByCoord.out.bam.bai) \
	    $(foreach sample,$(SAMPLES),bam/$(sample).bam) \
	    $(foreach sample,$(SAMPLES),bam/$(sample).bam.bai) \
	    summary/umi_summary.txt \
	    summary/umi_per_position_summary.txt

define merge-fastq
umi_tools/$1/$1_R1.fastq.gz : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G -w 72:00:00,"zcat $$(^) | gzip -c > $$(@)")
	
umi_tools/$1/$1_R2.fastq.gz : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G -w 72:00:00,"zcat $$(^) | gzip -c > $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merge-fastq,$(sample),$(split.$(sample)))))

define clip-fastq
umi_tools/$1/$1_R1_cl.fastq.gz : umi_tools/$1/$1_R1.fastq.gz umi_tools/$1/$1_R2.fastq.gz
	$$(call RUN,-c -n 1 -s 8G -m 16G -v $(UMITOOLS_ENV),"set -o pipefail && \
							     umi_tools extract \
							     -p $$(UMI_PATTERN) \
							     -I umi_tools/$1/$1_R1.fastq.gz \
							     -S umi_tools/$1/$1_R1_cl.fastq.gz \
							     --bc-pattern2=$$(UMI_PATTERN) \
							     --read2-in=umi_tools/$1/$1_R2.fastq.gz \
							     --read2-out=umi_tools/$1/$1_R2_cl.fastq.gz \
							     --log=umi_tools/$1/$1_cl.log")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call clip-fastq,$(sample))))

define align-fastq
star/$1.Aligned.sortedByCoord.out.bam : umi_tools/$1/$1_R1_cl.fastq.gz
	$$(call RUN,-n 4 -s 6G -m 10G -w 1440,"set -o pipefail && \
					       STAR $$(STAR_OPTS) \
					       --outFileNamePrefix star/$1. \
					       --runThreadN 4 \
					       --outSAMattrRGline \"ID:$1\" \"LB:$1\" \"SM:$1\" \"PL:$${SEQ_PLATFORM}\" \
					       --readFilesIn umi_tools/$1/$1_R1_cl.fastq.gz umi_tools/$1/$1_R2_cl.fastq.gz \
					       --readFilesCommand zcat")
                                   
star/$1.Aligned.sortedByCoord.out.bam.bai : star/$1.Aligned.sortedByCoord.out.bam
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
				      $$(SAMTOOLS) index $$(<)")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call align-fastq,$(sample))))

define dedup-bam
bam/$1.bam : star/$1.Aligned.sortedByCoord.out.bam star/$1.Aligned.sortedByCoord.out.bam.bai
	$$(call RUN,-c -n 1 -s 24G -m 48G -v $(UMITOOLS_ENV) -w 1440,"set -o pipefail && \
								      umi_tools dedup \
								      -I $$(<) \
								      -S $$(@) \
								      --output-stats=umi_tools/$1/$1 \
								      --paired")
																  
bam/$1.bam.bai : bam/$1.bam
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
				      $$(SAMTOOLS) index $$(<)")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call dedup-bam,$(sample))))
	
summary/umi_summary.txt : $(foreach sample,$(SAMPLES),bam/$(sample).bam)
	$(call RUN,-n 1 -s 48G -m 96G,"$(RSCRIPT) $(SCRIPTS_DIR)/rna_seq/summarize_umi.R --option 1 --sample_names '$(SAMPLES)'")
	
summary/umi_per_position_summary.txt : $(foreach sample,$(SAMPLES),bam/$(sample).bam)
	$(call RUN,-n 1 -s 48G -m 96G,"$(RSCRIPT) $(SCRIPTS_DIR)/rna_seq/summarize_umi.R --option 2 --sample_names '$(SAMPLES)'")
    
..DUMMY := $(shell mkdir -p version; \
	     $(UMITOOLS_ENV)/bin/umi_tools --version > version/umi_tools.txt; \
             echo "STAR" >> version/umi_tools.txt; \
	     STAR --version >> version/umi_tools.txt; \
             $(SAMTOOLS) --version >> version/umi_tools.txt)
.SECONDARY: 
.DELETE_ON_ERROR:
.PHONY: umi_tools
