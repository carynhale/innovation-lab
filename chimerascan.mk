# Chimerascan

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc

LOGDIR = log/chimscan.$(NOW)

CHIMSCAN_PYTHONPATH := /home/limr/share/usr/lib/python:/home/limr/share/usr/lib/python2.7
CHIMSCAN_PYTHON := $(HOME)/share/usr/bin/python
CHIMSCAN_INDEX := $(HOME)/share/reference/chimerascan_index
CHIMSCAN_NORMAL_FILTER = $(HOME)/share/scripts/normalFilterChimerascan.pl

RECURRENT_FUSIONS = $(RSCRIPT) $(HOME)/share/scripts/recurrentFusions.R

ONCOFUSE_MEM = $(JAVA) -Xmx$1 -jar $(HOME)/share/usr/oncofuse-v1.0.6/Oncofuse.jar
ONCOFUSE_TISSUE_TYPE ?= EPI

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all 

#ALL = $(foreach sample,$(SAMPLES),chimscan/$(sample).chimscan_timestamp)
ifdef NORMAL_CHIMSCAN_RESULTS
ALL += $(foreach sample,$(SAMPLES),chimscan/bedpe/$(sample).chimscan.nft.bedpe)
ALLTABLE = chimscan/alltables/all.chimscan.nft.oncofuse.merged.txt
ALL += chimscan/recur_tables/recurFusions.chimscan.nft.gene.txt
else 
ALLTABLE = chimscan/alltables/all.chimscan.oncofuse.merged.txt
ALL += chimscan/recur_tables/recurFusions.chimscan.nft.gene.txt
endif
ALL += $(ALLTABLE)

all : $(ALL)

CHIMERASCAN = PYTHONPATH=$(CHIMSCAN_PYTHONPATH) $(CHIMSCAN_PYTHON) /home/limr/share/usr/lib/python/chimerascan/chimerascan_run.py
CHIMERASCAN_OPTS = -v --quals illumina

chimscan/%.chimscan_timestamp : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call LSCRIPT_PARALLEL_MEM,4,6G,12G,"$(CHIMERASCAN) $(CHIMERASCAN_OPTS) -p 4 $(CHIMSCAN_INDEX) $^ $(@D)/$* && touch $@")

#chimerascan/tables/all.chimscan_results.txt : $(foreach sample,$(SAMPLES),chimerascan/$(sample).chimscan_timestamp)
#$(INIT) head -1 $(basename $<)/chimeras.bedpe > $@ && for x in $(addsuffix /chimeras.bedpe,$(basename $^)); do sed '1d' $$x >> $@; done

chimscan/bedpe/%.chimscan.bedpe : chimscan/%.chimscan_timestamp
	$(INIT) cp -f $(basename $<)/chimeras.bedpe $@ && rm -r chimscan/$*

%.chimscan.nft.bedpe : %.chimscan.bedpe
	$(call LSCRIPT_MEM,2G,4G,"$(PERL) $(CHIMSCAN_NORMAL_FILTER) -w 1000 $(NORMAL_CHIMSCAN_RESULTS) $< > $@")

chimscan/alltables/all.chimscan%txt : $(foreach sample,$(SAMPLES),chimscan/bedpe/$(sample).chimscan%bedpe)
	$(INIT) head -1 $< | sed 's/^/Sample\t/; s/#//' > $@ && for i in $^; do sed "1d; s/^/$$(basename $${i%%.chimscan_results.txt})\t/" $$i >> $@; done

chimscan/recur_tables/recurFusions.%.gene.txt : chimscan/alltables/all.%.txt
	$(INIT) $(RECURRENT_FUSIONS) --geneCol1 genes5p --geneCol2 genes3p --sampleCol Sample --outPrefix chimscan/recur_tables/recurGenes.$*  $< 

chimscan/alltables/all.chimscan%coord.txt : chimscan/alltables/all.chimscan%txt
	$(INIT) perl -lane 'if ($$. > 1) { $$coord5 = (($$F[9] eq "+")? $$F[3] + 1 : $$F[2] - 1) + 1; $$coord3 = (($$F[10] eq "+")? $$F[5] - 1 : $$F[6] + 1) + 1; print "$$F[1]\t$$coord5\t$$F[4]\t$$coord3\tEPI"; }' $< > $@

%.oncofuse.txt : %.coord.txt
	$(call LSCRIPT_MEM,8G,12G,"$(call ONCOFUSE_MEM,7G) $< coord $(ONCOFUSE_TISSUE_TYPE) $@")

%.oncofuse.merged.txt : %.txt %.oncofuse.txt 
	$(INIT) head -1 $< | sed 's/^/RowID\t/' > $<.tmp && awk 'BEGIN {OFS = "\t" } NR > 1 { print NR-1, $$0 }' $< >> $<.tmp ;\
		cut -f 2- $(<<) > $(<<).tmp; \
		$(RSCRIPT) $(MERGE) -X --byColX 1 --byColY 1 -H $<.tmp $(<<).tmp > $@ && rm -f $<.tmp $(<<).tmp

