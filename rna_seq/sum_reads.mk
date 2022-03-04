include innovation-lab/Makefile.inc

LOGDIR = log/sum_reads.$(NOW)

sum_reads : $(foreach sample,$(SAMPLES),sumreads/$(sample)_bygene.txt) \
	    sumreads/rpkm_bygene.txt \
	    sumreads/counts_bygene.txt
#	    $(foreach sample,$(SAMPLES),sumreads/$(sample)_byexon.txt)
#	    sumreads/rpkm_byexon.txt \
#           sumreads/counts_byexon.txt

define sum-reads
sumreads/$1_bygene.txt : bam/$1.bam bam/$1.bam.bai
	$$(call RUN,-n 1 -s 24G -m 48G -v $(SUMREADS_ENV),"$$(RSCRIPT) $$(SCRIPTS_DIR)/rna_seq/summarize_sumreads.R \
							   --option 1 \
							   --in_file $$(<) \
							   --out_file $$(@)")

#sumreads/$1_byexon.txt : bam/$1.bam bam/$1.bam.bai
#	$$(call RUN,-n 1 -s 24G -m 48G -v $(SUMREADS_ENV),"$$(RSCRIPT) $$(SCRIPTS_DIR)/rna_seq/summarize_sumreads.R \
#							   --option 2 \
#							   --in_file $$(<) \
#							   --out_file $$(@)")
							  
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call sum-reads,$(sample))))


sumreads/rpkm_bygene.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample)_bygene.txt)
	cut -f 1-2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 7 $$x | sed "s/_bygene/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

sumreads/counts_bygene.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample)_bygene.txt)
	cut -f 1-2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 3 $$x | sed "s/_bygene/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

#sumreads/rpkm_byexon.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads_byexon.txt)
#	cut -f 1-2 $< > $@; \
#	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 6 $$x | sed "s/exonRPKM/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done


#sumreads/counts_byexon.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads_byexon.txt)
#	cut -f 1-2 $< > $@; \
#	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 4 $$x | sed "s/exonCount/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done


..DUMMY := $(shell mkdir -p version; \
	     $(SUMREADS_ENV)/bin/R --version >> version/sum_reads.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: sum_reads
