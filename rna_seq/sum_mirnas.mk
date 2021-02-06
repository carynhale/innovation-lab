include innovation-lab/Makefile.inc

LOGDIR = log/sum_reads.$(NOW)

SUM_READS_RSCRIPT = $(SCRIPTS_DIR)/rna_seq/summarize_mirnas.R

sum_reads : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads_bymirna.txt) \
	    sumreads/rpkm_bymirna.txt \
	    sumreads/counts_bymirna.txt \
            

sumreads/%.sumreads_bymirna.txt : bam/%.bam bam/%.bam.bai
	$(call RUN,-n 1 -s 24G -m 48G -v $(SUMREADS_ENV),"$(RSCRIPT) $(SUM_READS_RSCRIPT) \
							  --in_file $(<) \
							  --out_file $(@)")

sumreads/rpkm_bymirna.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads_bymirna.txt)
	cut -f 1 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 4 $$x | sed "s/RPKM/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

sumreads/counts_bymirna.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads_bymirna.txt)
	cut -f 1 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 2 $$x | sed "s/ReadCount/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

..DUMMY := $(shell mkdir -p version; \
	     R --version >> version/sum_reads.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: sum_reads
