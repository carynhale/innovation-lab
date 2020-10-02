include innovation-lab/Makefile.inc

LOGDIR ?= log/defuse_umi.$(NOW)

defuse_umi : $(foreach sample,$(SAMPLES),defuse/$(sample).dedup.1.fastq) \
		 	 $(foreach sample,$(SAMPLES),defuse/$(sample).dedup/results.candidate.tsv) \
		 	 $(foreach sample,$(SAMPLES),defuse/$(sample).dedup/taskcomplete) \
		 	 defuse/summary.dedup.txt
		 
DEFUSE_CONFIG = innovation-lab/config/defuse.inc
DEFUSE_E75 = /home/brownd7/share/lib/resource_files/defuse/homo_sapiens/Ensembl/Grch37.p13/Sequence/defuse_e75

defuse/%.dedup.1.fastq : %.bam
	$$(call RUN,-n 4 -s 4G -m 9G,"set -o pipefail && \
								  $(SAMTOOLS) sort -T $(<D)/$* -O bam -n -@ 4 -m 6G $< | $(SAMTOOLS) fastq -f 1 -1 defuse/$*.dedup.1.fastq -2 defuse/$*.dedup.2.fastq")


define run-defuse
defuse/%.dedup/results.candidate.tsv : defuse/%.dedup.1.fastq defuse/%.dedup.2.fastq
	$$(call RUN,-c -n 10 -s 2G -m 3G -w 72:00:00 -v $(DEFUSE_ENV),"set -o pipefail && \
												  				   mkdir -p defuse && \
												  				   $$(DEFUSE) \
												  				   --config $$(DEFUSE_CONFIG) \
												  				   --dataset $$(DEFUSE_E75) \
												  				   --output defuse/$$(*).dedup \
												  				   --res defuse/$$(*).dedup/results.candidate.tsv \
												  				   --rescla defuse/$$(*).dedup/results.classify.tsv \
												  				   --resfil defuse/$$(*).dedup/results.filtered.tsv \
												  				   -1 defuse/$$(*).dedup.1.fastq \
												  				   -2 defuse/$$(*).dedup.2.fastq \
												  				   -s direct \
												  				   -p 10")
												  
defuse/%.dedup/taskcomplete : defuse/%.dedup/results.candidate.tsv
	$$(call RUN,-c -n 1 -s 1G -m 2G,"set -o pipefail && \
									 echo $$(*) > defuse/$$(*).dedup/taskcomplete")
	
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call run-defuse,$(sample))))
		
defuse/summary.dedup.txt : $(foreach sample,$(SAMPLES),defuse/$(sample)/taskcomplete)
	echo "cluster_id	splitr_sequence	splitr_count	splitr_span_pvalue	splitr_pos_pvalue	splitr_min_pvalue	adjacent	altsplice	break_adj_entropy1	break_adj_entropy2	break_adj_entropy_min	breakpoint_homology	breakseqs_estislands_percident	cdna_breakseqs_percident	deletion	est_breakseqs_percident	eversion	exonboundaries	expression1	expression2	gene1	gene2	gene_align_strand1	gene_align_strand2	gene_chromosome1	gene_chromosome2	gene_end1	gene_end2	gene_location1	gene_location2	gene_name1	gene_name2	gene_start1	gene_start2	gene_strand1	gene_strand2	genome_breakseqs_percident	genomic_break_pos1	genomic_break_pos2	genomic_ends1	genomic_ends2	genomic_starts1	genomic_starts2	genomic_strand1	genomic_strand2	interchromosomal	interrupted_index1	interrupted_index2	inversion	library_name	max_map_count	max_repeat_proportion	mean_map_count	min_map_count	num_multi_map	num_splice_variants	orf	read_through	repeat_proportion1	repeat_proportion2	span_count	span_coverage1	span_coverage2	span_coverage_max	span_coverage_min	splice_score	splicing_index1	splicing_index2	transcript1	transcript2" > defuse/summary.dedup.txt; \
	for i in $(SAMPLES); do \
		sed -e "1d" defuse/$$i.dedup/results.candidate.tsv | sed "s/$$//" >> defuse/summary.dedup.txt; \
	done
		
..DUMMY := $(shell mkdir -p version; \
			 ~/share/usr/env/defuse-0.8.0/bin/defuse_run.pl --help &> version/defuse_umi.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: defuse_umi
