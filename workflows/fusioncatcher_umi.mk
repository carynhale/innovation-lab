include innovation-lab/Makefile.inc

LOGDIR ?= log/fusioncatcher_umi.$(NOW)

fusioncatcher_umi : $(foreach sample,$(SAMPLES),fusioncatcher/$(sample).dedup/$(sample).1.fastq.gz)
#					$(foreach sample,$(SAMPLES),fusioncatcher/$(sample).dedup/out/taskcomplete) \
#					fusioncatcher/summary.dedup.txt

CACHE = $(HOME)/share/usr/env/fusioncatcher-1.2.0/share/fusioncatcher-1.20/db/current

define fusioncatcher-dedup
fusioncatcher/$1.dedup/$1.1.fastq.gz : bam/$1.bam
	$$(call RUN,-n 4 -s 4G -m 9G,"set -o pipefail && \
								  $$(SAMTOOLS) sort -T $$(<D)/$$(*) -O bam -n -@ 4 -m 6G $$(<) | \
								  $$(SAMTOOLS) fastq -f 1 -1 fusioncatcher/$$(*).dedup/$$(*).1.fastq -2 fusioncatcher/$$(*).dedup/$$(*).2.fastq && \
								  gzip fusioncatcher/$$(*).dedup/$$(*).1.fastq && \
								  gzip fusioncatcher/$$(*).dedup/$$(*).2.fastq")

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
		
fusioncatcher/summary.dedup.txt : $(foreach sample,$(SAMPLES),fusioncatcher/$(sample).dedup/out/taskcomplete)
	echo "Gene_1_symbol(5end_fusion_partner)	Gene_2_symbol(3end_fusion_partner)	Fusion_description	Counts_of_common_mapping_reads	Spanning_pairs	Spanning_unique_reads	Longest_anchor_found	Fusion_finding_method	Fusion_point_for_gene_1(5end_fusion_partner)	Fusion_point_for_gene_2(3end_fusion_partner)	Gene_1_id(5end_fusion_partner)	Gene_2_id(3end_fusion_partner)	Exon_1_id(5end_fusion_partner)	Exon_2_id(3end_fusion_partner)	Fusion_sequence	Predicted_effect	Sample_name" > fusioncatcher/summary.dedup.txt; \
	for i in $(SAMPLES); do \
		sed -e "1d" fusioncatcher/$$i.dedup/out/final-list_candidate-fusion-genes.hg19.txt | sed "s/$$/\t$$i/" >> fusioncatcher/summary.dedup.txt; \
	done

..DUMMY := $(shell mkdir -p version; \
			 ~/share/usr/env/fusioncatcher-1.2.0/bin/fusioncatcher.py --version &> version/fusioncatcher_umi.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: fusioncatcher_umi
