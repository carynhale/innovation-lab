include innovation-lab/Makefile.inc
include innovation-lab/config/align.inc
include innovation-lab/config/gatk.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/fusioncatcher_umi.$(NOW)

fusioncatcher_umi : $(foreach sample,$(SAMPLES),umi_tools/$(sample)/$(sample)_R1.fastq.gz) \
					$(foreach sample,$(SAMPLES),umi_tools/$(sample)/$(sample)_R1_cl.fastq.gz) \
					$(foreach sample,$(SAMPLES),star/$(sample).Aligned.sortedByCoord.out.bam) \
					$(foreach sample,$(SAMPLES),star/$(sample).Aligned.sortedByCoord.out.bam.bai) \
					$(foreach sample,$(SAMPLES),bam/$(sample).bam) \
					$(foreach sample,$(SAMPLES),bam/$(sample).bam.bai) \
					$(foreach sample,$(SAMPLES),fusioncatcher/$(sample)/$(sample).1.fastq.gz) \
					$(foreach sample,$(SAMPLES),fusioncatcher/$(sample)/out/taskcomplete) \
					fusioncatcher/summary.txt \
					$(foreach sample,$(SAMPLES),fusioncatcher/$(sample).dedup/$(sample).1.fastq.gz) \
					$(foreach sample,$(SAMPLES),fusioncatcher/$(sample).dedup/out/taskcomplete) \
					fusioncatcher/summary.dedup.txt

#ALIGNER := star
#BAM_NO_REALN = true
#BAM_NO_RECAL = true
#BAM_NO_FILTER = true
#BAM_DUP_TYPE = none
#SEQ_PLATFROM = illumina
#
#UMI_PATTERN = NNNX
#
#STAR_OPTS = --genomeDir $(STAR_REF) \
#            --outSAMtype BAM SortedByCoordinate \
#			--twopassMode Basic \
#            --outReadsUnmapped None \
#            --outSAMunmapped Within \
#            --chimSegmentMin 12 \
#            --chimJunctionOverhangMin 12 \
#            --alignSJDBoverhangMin 10 \
#            --alignMatesGapMax 200000 \
#            --alignIntronMax 200000 \
#            --chimSegmentReadGapMax parameter 3 \
#            --alignSJstitchMismatchNmax 5 -1 5 5 \
#            --chimOutType WithinBAM \
#			--quantMode GeneCounts

#CACHE = $(HOME)/share/usr/env/fusioncatcher-1.2.0/share/fusioncatcher-1.20/db/current

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
	
define fusioncatcher
fusioncatcher/$1/$1.1.fastq.gz : star/$1.Aligned.sortedByCoord.out.bam
	$$(call RUN,-n 4 -s 4G -m 9G,"set -o pipefail && \
								  $$(SAMTOOLS) sort -T star/$1 -O bam -n -@ 4 -m 6G $$(<) | \
								  bedtools bamtofastq -i - -fq >(gzip -c > fusioncatcher/$1/$1.1.fastq.gz) -fq2 >(gzip -c > fusioncatcher/$1/$1.2.fastq.gz)")

fusioncatcher/$1/out/taskcomplete : fusioncatcher/$1/$1.1.fastq.gz
	$$(call RUN,-c -n 8 -s 2G -m 3G -v $(FUSIONCATCHER_ENV) -w 72:00:00,"set -o pipefail && \
																		 mkdir -p fusioncatcher/$1/out && \
																		 fusioncatcher.py \
																		 -i fusioncatcher/$1 \
																		 -o fusioncatcher/$1/out \
																		 -d $$(CACHE) \
																		 -p 8 && \
																		 echo $1 > fusioncatcher/$1/out/taskcomplete")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fusioncatcher,$(sample))))


define fusioncatcher-dedup
fusioncatcher/$1.dedup/$1.1.fastq.gz : bam/$1.bam
	$$(call RUN,-n 4 -s 4G -m 9G,"set -o pipefail && \
								  $$(SAMTOOLS) sort -T bam/$1 -O bam -n -@ 4 -m 6G $$(<) | \
								  bedtools bamtofastq -i - -fq >(gzip -c > fusioncatcher/$1.dedup/$1.1.fastq.gz) -fq2 >(gzip -c > fusioncatcher/$1.dedup/$1.2.fastq.gz)")

fusioncatcher/$1.dedup/out/taskcomplete : fusioncatcher/$1.dedup/$1.1.fastq.gz
	$$(call RUN,-c -n 8 -s 2G -m 3G -v $(FUSIONCATCHER_ENV) -w 72:00:00,"set -o pipefail && \
																		 mkdir -p fusioncatcher/$1.dedup/out && \
																		 fusioncatcher.py \
																		 -i fusioncatcher/$1.dedup \
																		 -o fusioncatcher/$1.dedup/out \
																		 -d $$(CACHE) \
																		 -p 8 && \
																		 echo $1 > fusioncatcher/$1.dedup/out/taskcomplete")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fusioncatcher-dedup,$(sample))))


fusioncatcher/summary.txt : $(foreach sample,$(SAMPLES),fusioncatcher/$(sample)/out/taskcomplete)
	echo "Gene_1_symbol(5end_fusion_partner)	Gene_2_symbol(3end_fusion_partner)	Fusion_description	Counts_of_common_mapping_reads	Spanning_pairs	Spanning_unique_reads	Longest_anchor_found	Fusion_finding_method	Fusion_point_for_gene_1(5end_fusion_partner)	Fusion_point_for_gene_2(3end_fusion_partner)	Gene_1_id(5end_fusion_partner)	Gene_2_id(3end_fusion_partner)	Exon_1_id(5end_fusion_partner)	Exon_2_id(3end_fusion_partner)	Fusion_sequence	Predicted_effect	Sample_name" > fusioncatcher/summary.txt; \
	for i in $(SAMPLES); do \
		sed -e "1d" fusioncatcher/$$i/out/final-list_candidate-fusion-genes.hg19.txt | sed "s/$$/\t$$i/" >> fusioncatcher/summary.txt; \
	done


fusioncatcher/summary.dedup.txt : $(foreach sample,$(SAMPLES),fusioncatcher/$(sample).dedup/out/taskcomplete)
	echo "Gene_1_symbol(5end_fusion_partner)	Gene_2_symbol(3end_fusion_partner)	Fusion_description	Counts_of_common_mapping_reads	Spanning_pairs	Spanning_unique_reads	Longest_anchor_found	Fusion_finding_method	Fusion_point_for_gene_1(5end_fusion_partner)	Fusion_point_for_gene_2(3end_fusion_partner)	Gene_1_id(5end_fusion_partner)	Gene_2_id(3end_fusion_partner)	Exon_1_id(5end_fusion_partner)	Exon_2_id(3end_fusion_partner)	Fusion_sequence	Predicted_effect	Sample_name" > fusioncatcher/summary.dedup.txt; \
	for i in $(SAMPLES); do \
		sed -e "1d" fusioncatcher/$$i.dedup/out/final-list_candidate-fusion-genes.hg19.txt | sed "s/$$/\t$$i/" >> fusioncatcher/summary.dedup.txt; \
	done


..DUMMY := $(shell mkdir -p version; \
			 ~/share/usr/env/fusioncatcher-1.2.0/bin/fusioncatcher.py --version &> version/fusioncatcher_umi.txt; \
			 $(UMITOOLS_ENV)/bin/umi_tools --version >> version/fusioncatcher_umi.txt; \
             echo "STAR" >> version/fusioncatcher_umi.txt; \
			 STAR --version >> version/fusioncatcher_umi.txt; \
             $(SAMTOOLS) --version >> version/fusioncatcher_umi.txt)
.SECONDARY: 
.DELETE_ON_ERROR:
.PHONY: fusioncatcher_umi

#include innovation-lab/fastq_tools/fastq.mk
#include innovation-lab/bam_tools/process_bam.mk
#include innovation-lab/aligners/align.mk