include innovation-lab/Makefile.inc

LOGDIR ?= log/defuse.$(NOW)

defuse : $(foreach sample,$(SAMPLES),defuse/$(sample)/$(sample).1.fastq) \
	 $(foreach sample,$(SAMPLES),defuse/$(sample)/$(sample).2.fastq) \
	 $(foreach sample,$(SAMPLES),defuse/$(sample)/results.candidate.tsv) \
	 $(foreach sample,$(SAMPLES),defuse/$(sample)/taskcomplete) \
	 defuse/summary.txt
		 
DEFUSE_CONFIG = innovation-lab/config/defuse.inc
DEFUSE_E75 = /home/brownd7/share/lib/resource_files/defuse/homo_sapiens/Ensembl/Grch37.p13/Sequence/defuse_e75

define merged-fastq
defuse/$1/$1.1.fastq : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 12 -s 0.5G -m 1G -v $(PIGZ_ENV),"set -o pipefail && \
							$$(PIGZ) -cd -p 12 $$(^) > $$(@)")
					 
defuse/$1/$1.2.fastq : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 12 -s 0.5G -m 1G -v $(PIGZ_ENV),"set -o pipefail && \
							$$(PIGZ) -cd -p 12 $$(^) > $$(@)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merged-fastq,$(sample),$(split.$(sample)))))

define run-defuse
defuse/$1/results.candidate.tsv : defuse/$1/$1.1.fastq defuse/$1/$1.2.fastq
	$$(call RUN,-c -n 10 -s 2G -m 3G -w 72:00:00 -v $(DEFUSE_ENV),"set -o pipefail && \
								       mkdir -p defuse && \
								       $$(DEFUSE) \
								       --config $$(DEFUSE_CONFIG) \
								       --dataset $$(DEFUSE_E75) \
								       --output defuse/$$(*) \
								       --res defuse/$$(*)/results.candidate.tsv \
								       --rescla defuse/$$(*)/results.classify.tsv \
								       --resfil defuse/$$(*)/results.filtered.tsv \
								       -1 $$(<) \
								       -2 $$(<<) \
								       -s direct \
								       -p 10")
												  
defuse/$1/taskcomplete : defuse/$1/results.candidate.tsv
	$$(call RUN,-c -n 1 -s 1G -m 2G,"set -o pipefail && \
					 echo $$(*) > $$(@)")
	
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call run-defuse,$(sample))))
		
defuse/summary.txt : $(foreach sample,$(SAMPLES),defuse/$(sample)/taskcomplete)
	echo "cluster_id	splitr_sequence	splitr_count	splitr_span_pvalue	splitr_pos_pvalue	splitr_min_pvalue	adjacent	altsplice	break_adj_entropy1	break_adj_entropy2	break_adj_entropy_min	breakpoint_homology	breakseqs_estislands_percident	cdna_breakseqs_percident	deletion	est_breakseqs_percident	eversion	exonboundaries	expression1	expression2	gene1	gene2	gene_align_strand1	gene_align_strand2	gene_chromosome1	gene_chromosome2	gene_end1	gene_end2	gene_location1	gene_location2	gene_name1	gene_name2	gene_start1	gene_start2	gene_strand1	gene_strand2	genome_breakseqs_percident	genomic_break_pos1	genomic_break_pos2	genomic_ends1	genomic_ends2	genomic_starts1	genomic_starts2	genomic_strand1	genomic_strand2	interchromosomal	interrupted_index1	interrupted_index2	inversion	library_name	max_map_count	max_repeat_proportion	mean_map_count	min_map_count	num_multi_map	num_splice_variants	orf	read_through	repeat_proportion1	repeat_proportion2	span_count	span_coverage1	span_coverage2	span_coverage_max	span_coverage_min	splice_score	splicing_index1	splicing_index2	transcript1	transcript2" > defuse/summary.txt; \
	for i in $(SAMPLES); do \
		sed -e "1d" defuse/$$i/results.candidate.tsv | sed "s/$$//" >> defuse/summary.txt; \
	done
		
..DUMMY := $(shell mkdir -p version; \
	     ~/share/usr/env/defuse-0.8.0/bin/defuse_run.pl --help &> version/defuse.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: defuse
