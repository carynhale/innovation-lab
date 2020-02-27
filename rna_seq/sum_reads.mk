include innovation-lab/Makefile.inc

LOGDIR = log/sum_reads.$(NOW)

SUM_READS_RSCRIPT = $(SCRIPTS_DIR)/rna_seq/summarize_rnaseqreads.R
SUM_EXONS_RSCRIPT = $(SCRIPTS_DIR)/rna_seq/summarize_rnaseqreadsbyexon.R

sum_reads : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads_bygene.txt) \
            $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads_byexon.txt) \
            sumreads/rpkm_bygene.txt \
            sumreads/rpkm_byexon.txt \
            sumreads/counts_bygene.txt \
            sumreads/counts_byexon.txt

sumreads/%.sumreads_bygene.txt : bam/%.bam bam/%.bam.bai
	$(call RUN,-n 1 -s 24G -m 48G -v $(SUMREADS_ENV),"$(RSCRIPT) $(SUM_READS_RSCRIPT) \
                                                      --genome $(REF) \
                                                      --in_file $(<) \
                                                      --out_file $(@)")

sumreads/%.sumreads_byexon.txt : bam/%.bam bam/%.bam.bai
	$(call RUN,-n 1 -s 24G -m 48G -v $(SUMREADS_ENV),"$(RSCRIPT) $(SUM_EXONS_RSCRIPT) \
                                                      --genome $(REF) \
                                                      --in_file $(<) \
                                                      --out_file $(@)")

sumreads/rpkm_bygene.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads_bygene.txt)
	cut -f 2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 7 $$x | sed "s/exonRPKM/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

sumreads/rpkm_byexon.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads_byexon.txt)
	cut -f 1-2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 6 $$x | sed "s/exonRPKM/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

sumreads/counts_bygene.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads_bygene.txt)
	cut -f 2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 3 $$x | sed "s/countsByGene/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done

sumreads/counts_byexon.txt : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads_byexon.txt)
	cut -f 1-2 $< > $@; \
	for x in $^; do sample=`echo $$x | sed 's/.*\///; s/\..*//'`; cut -f 4 $$x | sed "s/exonCount/$$sample/" | paste $@ - > $@.tmp; mv $@.tmp $@; done


..DUMMY := $(shell mkdir -p version; \
			 R --version >> version/sum_reads.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: sum_reads
