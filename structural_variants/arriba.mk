include innovation-lab/Makefile.inc
include innovation-lab/config/arriba.inc

LOGDIR ?= log/arriba.$(NOW)

arriba : $(foreach sample,$(SAMPLES),arriba/$(sample)/$(sample).1.fastq.gz) \
	 $(foreach sample,$(SAMPLES),arriba/$(sample)/$(sample).2.fastq.gz) \
	 $(foreach sample,$(SAMPLES),arriba/$(sample)/fusions.tsv) \
	 $(foreach sample,$(SAMPLES),arriba/$(sample)/fusions.pdf) \
	 arriba/summary.txt
	 
	 
define merged-fastq
arriba/$1/$1.1.fastq.gz : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 12 -s 0.5G -m 1G -v $(PIGZ_ENV),"set -o pipefail && \
							$$(PIGZ) -cd -p 12 $$(^) | $$(PIGZ) -c -p 12 > $$(@)")
					 
arriba/$1/$1.2.fastq.gz : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 12 -s 0.5G -m 1G -v $(PIGZ_ENV),"set -o pipefail && \
							$$(PIGZ) -cd -p 12 $$(^) | $$(PIGZ) -c -p 12 > $$(@)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merged-fastq,$(sample),$(split.$(sample)))))

define arriba
arriba/$1/$1.Aligned.out.bam : arriba/$1/$1.1.fastq.gz arriba/$1/$1.2.fastq.gz
	$$(call RUN,-c -n $(STAR_THREADS) -s 3G -m 4G -v $(ARRIBA_ENV) -w 72:00:00,"set -o pipefail && \
										    STAR \
										    --runThreadN $$(STAR_THREADS) \
										    --genomeDir $$(STAR_INDEX_DIR) \
										    --genomeLoad NoSharedMemory \
										    --readFilesIn $$(<) $$(<<) \
										    --readFilesCommand zcat \
										    --outStd BAM_Unsorted \
										    --outSAMtype BAM Unsorted \
										    --outSAMunmapped Within \
										    --outBAMcompression 0 \
										    --outFilterMultimapNmax 50 \
										    --peOverlapNbasesMin 10 \
										    --alignSplicedMateMapLminOverLmate 0.5 \
										    --alignSJstitchMismatchNmax 5 -1 5 5 \
										    --chimSegmentMin 10 \
										    --chimOutType WithinBAM HardClip \
										    --chimJunctionOverhangMin 10 \
										    --chimScoreDropMax 30 \
										    --chimScoreJunctionNonGTAG 0 \
										    --chimScoreSeparation 1 \
										    --chimSegmentReadGapMax 3 \
										    --chimMultimapNmax 50 \
										    --outFileNamePrefix arriba/$1/$1. > arriba/$1/$1.Aligned.out.bam")


arriba/$1/fusions.tsv : arriba/$1/$1.Aligned.out.bam
	$$(call RUN,-c -n 1 -s 24G -m 36G -v $(ARRIBA_ENV),"set -o pipefail && \
							    $$(ARRIBA_EXE) -x arriba/$1/$1.Aligned.out.bam \
							    -o arriba/$1/fusions.tsv \
							    -O arriba/$1/discarded.tsv \
							    -a $$(ASSEMBLY_FA) \
							    -g $$(ANNOTATION_GTF) \
							    -b $$(BLACKLIST_TSV) \
							    -k $$(KNOWN_FUSIONS_TSV) \
							    -t $$(KNOWN_FUSIONS_TSV) \
							    -p $$(PROTEIN_DOMAINS_GFF3)")
							    
arriba/$1/fusions.pdf : arriba/$1/fusions.tsv arriba/$1/$1.Aligned.out.bam
	$$(call RUN,-c -n 1 -s 12G -m 24G -v $(GENOMIC_ALIGNMENTS_ENV),"set -o pipefail && \
									$$(RSCRIPT) $$(DRAW_FUSIONS) \
									--fusions=$$(<) \
									--annotation=$$(ANNOTATION_GTF) \
									--alignments=$$(<<) \
									--cytobands=$$(CYTOBAND) \
									--proteinDomains=$$(PROTEIN_DOMAINS_GFF3) \
									--output=$$(@)")


endef
$(foreach sample,$(SAMPLES),\
                $(eval $(call arriba,$(sample))))
				
arriba/summary.txt : $(foreach sample,$(SAMPLES),arriba/$(sample)/fusions.tsv)
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
					   $(RSCRIPT) $(SCRIPTS_DIR)/rna_seq/summarize_arriba.R --sample_names '$(SAMPLES)'")
									   
..DUMMY := $(shell mkdir -p version; \
	     $(ARRIBA_EXE) -h > version/arriba.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: arriba
