include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/bam_interval_metrics.$(NOW)
PHONY += metrics metrics/pileup metrics/cov metrics/picard metrics/summary metrics/report

TARGETS_LIST ?= $(HOME)/share/reference/target_panels/MSK-IMPACT-468.list

interval_metrics : $(foreach sample,$(SAMPLES),metrics/pileup/$(sample)-pileup.txt) \
			   	   $(foreach sample,$(SAMPLES),metrics/cov/$(sample)-ontarget.txt) \
				   $(foreach sample,$(SAMPLES),metrics/cov/$(sample)-offtarget.txt) \
				   $(foreach sample,$(SAMPLES),metrics/picard/$(sample)-idx_stats.txt) \
				   $(foreach sample,$(SAMPLES),metrics/picard/$(sample)-aln_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/picard/$(sample)-insert_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/picard/$(sample)-oxog_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/picard/$(sample)-hs_metrics.txt) \
				   metrics/summary/metrics_idx.tsv \
				   metrics/summary/metrics_aln.tsv \
				   metrics/summary/metrics_insert.tsv \
				   metrics/summary/metrics_insert_distribution.tsv \
				   metrics/summary/metrics_oxog.tsv \
				   metrics/summary/metrics_hs.tsv \
				   metrics/summary/metrics_coverage.tsv \
				   metrics/report/target_coverage.pdf \
				   metrics/report/alignment_summary.pdf \
				   metrics/report/insert_size_summary.pdf \
				   metrics/report/insert_size_distribution.pdf \
				   metrics/report/non_reference_calls.pdf
				   
define pileup-metric
metrics/pileup/$1-pileup.txt : bam/$1.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
								 cp bam/$1.bam metrics/pileup/$1.bam && \
								 cp bam/$1.bam.bai metrics/pileup/$1.bam.bai && \
								 cp bam/$1.bai metrics/pileup/$1.bai && \
								 cd metrics/pileup && \
								 $(HOME)/share/usr/jdk1.8.0_74/bin/java -server -Xms2G -Xmx8G -cp $(WALTZ) org.mskcc.juber.waltz.Waltz PileupMetrics 20 $1.bam $(REF_FASTA) $(TARGETS_FILE) && \
								 rm -rf $1.bam && \
								 rm -rf $1.bam.bai && \
								 rm -rf $1.bai && \
								 cd ../..")
									 
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call pileup-metric,$(sample))))
				   
define coverage-metric
metrics/cov/$1-ontarget.txt : bam/$1.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
								 samtools view -L $$(ONTARGET_FILE) $$(<) -b > metrics/cov/$1-ontarget.bam && \
								 samtools index metrics/cov/$1-ontarget.bam && \
								 java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) BamIndexStats \
								 I=metrics/cov/$1-ontarget.bam \
								 TMP_DIR=$(TMPDIR) \
								 > $$(@) && \
								 rm -rf metrics/cov/$1-ontarget.bam && \
								 rm -rf metrics/cov/$1-ontarget.bam.bai")
									 
metrics/cov/$1-offtarget.txt : bam/$1.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
								 samtools view -L $$(OFFTARGET_FILE) $$(<) -b > metrics/cov/$1-offtarget.bam && \
								 samtools index metrics/cov/$1-offtarget.bam && \
								 java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) BamIndexStats \
								 I=metrics/cov/$1-offtarget.bam \
								 TMP_DIR=$(TMPDIR) \
								 > $$(@) && \
								 rm -rf metrics/cov/$1-offtarget.bam && \
								 rm -rf metrics/cov/$1-offtarget.bam.bai")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call coverage-metric,$(sample))))

define picard-metric
metrics/picard/$1-idx_stats.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) BamIndexStats \
									   I=$$(<) \
									   TMP_DIR=$(TMPDIR) \
									   > $$(@)")
									   
metrics/picard/$1-aln_metrics.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectAlignmentSummaryMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   TMP_DIR=$(TMPDIR)")

metrics/picard/$1-insert_metrics.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectInsertSizeMetrics \
									   I=$$(<) \
									   O=$$(@) \
									   H=metrics/picard/$1.insert_metrics.pdf \
									   M=0.5 \
									   TMP_DIR=$(TMPDIR)")
												
metrics/picard/$1-oxog_metrics.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectOxoGMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   TMP_DIR=$(TMPDIR)")

metrics/picard/$1-hs_metrics.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CalculateHsMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   BAIT_INTERVALS=$(TARGETS_LIST) \
									   TARGET_INTERVALS=$(TARGETS_LIST) \
									   TMP_DIR=$(TMPDIR)")
												
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call picard-metric,$(sample))))
		
		
metrics/summary/metrics_idx.tsv : $(wildcard metrics/picard/$(SAMPLES)-idx_stats.txt)
	$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/qc/bam_interval_metrics_summary.R --metric 1 --samples '$(SAMPLES)'")
									  
metrics/summary/metrics_aln.tsv : $(wildcard metrics/picard/$(SAMPLES)-aln_metrics.txt)
	$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/qc/bam_interval_metrics_summary.R --metric 2 --samples '$(SAMPLES)'")
									  
metrics/summary/metrics_insert.tsv : $(wildcard metrics/picard/$(SAMPLES)-insert_metrics.txt)
	$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/qc/bam_interval_metrics_summary.R --metric 3 --samples '$(SAMPLES)'")
									  
metrics/summary/metrics_insert_distribution.tsv : $(wildcard metrics/picard$(SAMPLES)-insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/qc/bam_interval_metrics_summary.R --metric 4 --samples '$(SAMPLES)'")
									  
metrics/summary/metrics_oxog.tsv : $(wildcard metrics/picard/$(SAMPLES)-oxog_metrics.txt)
	$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/qc/bam_interval_metrics_summary.R --metric 5 --samples '$(SAMPLES)'")									  
		
metrics/summary/metrics_hs.tsv : $(wildcard metrics/picard/$(SAMPLES)-hs_metrics.txt)
	$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/qc/bam_interval_metrics_summary.R --metric 6 --samples '$(SAMPLES)'")
									  
metrics/summary/metrics_coverage.tsv : $(wildcard metrics/cov/$(SAMPLES)-ontarget.txt) $(wildcard metrics/cov/$(SAMPLES)-offtarget.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/qc/bam_interval_metrics_summary.R --metric 7 --samples '$(SAMPLES)'")
									  
metrics/report/target_coverage.pdf : metrics/summary/metrics_hs.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/qc/bam_interval_metrics_plot.R --type 1")

metrics/report/alignment_summary.pdf : metrics/summary/metrics_aln.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/qc/bam_interval_metrics_plot.R --type 2")

metrics/report/insert_size_summary.pdf : metrics/summary/metrics_insert.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/qc/bam_interval_metrics_plot.R --type 3")
									  
metrics/report/insert_size_distribution.pdf : metrics/summary/metrics_insert_distribution.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/qc/bam_interval_metrics_plot.R --type 4")
									  
metrics/report/non_reference_calls.pdf : $(wildcard metrics/pileup/$(SAMPLES)-pileup.txt)
	$(call RUN, -c -n 1 -s 48G -m 72G -v $(SUPERHEAT_ENV),"set -o pipefail && \
														   $(RSCRIPT) modules/qc/bam_interval_metrics_plot.R --type 5 --sample_names '$(SAMPLES)' && \
														   gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dFirstPage=2 -dLastPage=2 -sOutputFile=metrics/report/non_reference_calls-2.pdf metrics/report/non_reference_calls.pdf && \
														   rm metrics/report/non_reference_calls.pdf && \
														   mv metrics/report/non_reference_calls-2.pdf metrics/report/non_reference_calls.pdf")
									  
									  		
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
