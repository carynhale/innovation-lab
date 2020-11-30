include innovation-lab/Makefile.inc

LOGDIR ?= log/arriba.$(NOW)

arriba : $(foreach sample,$(SAMPLES),arriba/$(sample)/fusion.tsv)

ARRIBA_ENV = $(HOME)/share/usr/env/arriba-2.0.0
ARRIBA_EXE = $(HOME)/share/usr/env/arriba-2.0.0/src/arriba_v2.0.0/arriba
STAR_INDEX_DIR = $(HOME)/share/usr/env/arriba-2.0.0/src/arriba_v2.0.0/STAR_index_GRCh37viral_GENCODE19
ANNOTATION_GTF = $(HOME)/share/usr/env/arriba-2.0.0/src/arriba_v2.0.0/GENCODE19.gtf
ASSEMBLY_FA = $(HOME)/share/usr/env/arriba-2.0.0/src/arriba_v2.0.0/GRCh37viral.fa
BLACKLIST_TSV = $(HOME)/share/usr/env/arriba-2.0.0/src/arriba_v2.0.0/database/blacklist_hg19_hs37d5_GRCh37_v2.0.0.tsv.gz
KNOWN_FUSIONS_TSV = $(HOME)/share/usr/env/arriba-2.0.0/src/arriba_v2.0.0/database/known_fusions_hg19_hs37d5_GRCh37_v2.0.0.tsv.gz
PROTEIN_DOMAINS_GFF3 = $(HOME)/share/usr/env/arriba-2.0.0/src/arriba_v2.0.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.0.0.gff3
THREADS = 16

define run-arriba
arriba/$1-nodedup/$1.Aligned.out.bam : umi_tools/$1/$1_R1_cl.fastq.gz
	$$(call RUN,-c -n $(THREADS) -s 1G -m 2G -v $(ARRIBA_ENV) -w 72:00:00,"set -o pipefail && \
										STAR \
										--runThreadN $$(THREADS) \
										--genomeDir $$(STAR_INDEX_DIR) \
										--genomeLoad NoSharedMemory \
										--readFilesIn umi_tools/$1/$1_R1_cl.fastq.gz umi_tools/$1/$1_R2_cl.fastq.gz \
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
										--outFileNamePrefix arriba/$1-nodedup/$1. > arriba/$1-nodedup/$1.Aligned.out.bam")


arriba/$1-nodedup/$1-nodedup.tsv : arriba/$1-nodedup/$1.Aligned.out.bam
	$$(call RUN,-c -n 1 -s 24G -m 36G -v $(ARRIBA_ENV),"set -o pipefail && \
								$$(ARRIBA_EXE) -x arriba/$1-nodedup/$1.Aligned.out.bam \
								-o arriba/$1-nodedup/$1-nodedup.tsv -O arriba/$1-nodedup/$1-nodedup.discarded.tsv \
								-a $$(ASSEMBLY_FA) \
								-g $$(ANNOTATION_GTF) \
								-b $$(BLACKLIST_TSV) \
								-k $$(KNOWN_FUSIONS_TSV) \
								-t $$(KNOWN_FUSIONS_TSV) \
								-p $$(PROTEIN_DOMAINS_GFF3)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call run-arriba,$(sample))))

..DUMMY := $(shell mkdir -p version; \
			 $(ARRIBA_EXE) -h > version/arriba.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: arriba
