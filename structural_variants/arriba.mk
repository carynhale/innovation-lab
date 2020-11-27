include innovation-lab/Makefile.inc

LOGDIR ?= log/arriba.$(NOW)

arriba : $(foreach sample,$(SAMPLES),arriba/$(sample)/$(sample).fusion.tsv) \
		 arriba/summary.txt
		 
DBSNP_SUBSET = $(HOME)/share/lib/bed_files/dbsnp_137.b37_subset.bed

STAR_INDEX_DIR = $(STAR_REF)
ASSEMBLY_FA = $(REF_FASTA)
ANNOTATION_GTF = $(HOME)/share/lib/resource_files/GENCODE19.gtf
BLACKLIST_TSV = $(HOME)/share/usr/arriba_v2.0.0/database/blacklist_hg19_hs37d5_GRCh37_v2.0.0.tsv
KNOWN_FUSIONS_TSV = $(HOME)/share/usr/arriba_v2.0.0/database/known_fusions_hg19_hs37d5_GRCh37_v2.0.0.tsv
PROTEIN_DOMAINS_GFF3 = $(HOME)/share/usr/arriba_v2.0.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.0.0.gff3
ARRIBA = $(HOME)/share/usr/arriba_v2.0.0/arriba

define run-arriba
arriba/$1/$1.fusion.tsv : bam/$1.bam
	$$(call RUN,-c -n 1 -s 8G -m 12G -v $(STARFUSION_ENV) -w 36:00:00,"set -o pipefail && \
																	  tee bam/$1.bam | \
																	  /home/brownd7/share/usr/arriba_v2.0.0/arriba \
																	  -x /dev/stdin \
																	  -o arriba/$1/$1.fusion.tsv -O arriba/$1/$1.discarded.tsv \
																	  -a $$(ASSEMBLY_FA) \
																	  -g $$(ANNOTATION_GTF) \
																	  -b $$(BLACKLIST_TSV) \
																	  -k $$(KNOWN_FUSIONS_TSV) \
																	  -t $$(KNOWN_FUSIONS_TSV) \
																	  -p $$(PROTEIN_DOMAINS_GFF3)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call run-arriba,$(sample))))
		
arriba/summary.txt : $(foreach sample,$(SAMPLES),arriba/$(sample)/$(sample).fusion.tsv)
	echo "" > arriba/summary.txt; \
	for i in $(SAMPLES); do \
		sed -e "1d" arriba/$$i/$$i.fusion.tsv | sed "s/$$/\t$$i/" >> arriba/summary.txt; \
	done

..DUMMY := $(shell mkdir -p version; \
			 $(ARRIBA) -h > version/arriba.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: arriba
