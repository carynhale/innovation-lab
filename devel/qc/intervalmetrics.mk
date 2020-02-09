include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/interval_metrics.$(NOW)
PHONY += metrics metrics/standard metrics/unfiltered metrics/duplex metrics/simplex metrics/summary

POOL_A_TARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.list
POOL_B_TARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.list
ONTARGET_FILE_A ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.bed
ONTARGET_FILE_B ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.bed
OFFTARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-AB.offtarget.bed

interval_metrics : $(foreach sample,$(SAMPLES),metrics/standard/$(sample).idx_stats.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).aln_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).insert_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).oxog_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).probe-A.hs_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).probe-B.hs_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).probe-A.hs_metrics-nodedup.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).probe-B.hs_metrics-nodedup.txt) \
				   metrics/standard/metrics_idx.tsv \
				   metrics/standard/metrics_aln.tsv \
				   metrics/standard/metrics_insert.tsv \
				   metrics/standard/metrics_insert_distribution.tsv \
				   metrics/standard/metrics_oxog.tsv \
				   metrics/standard/metrics_hs.tsv \
				   $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample).idx_stats.txt) \
				   $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample).aln_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample).insert_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample).probe-A.hs_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample).probe-B.hs_metrics.txt) \
				   metrics/unfiltered/metrics_idx.tsv \
				   metrics/unfiltered/metrics_aln.tsv \
				   metrics/unfiltered/metrics_insert.tsv \
				   metrics/unfiltered/metrics_insert_distribution.tsv \
				   metrics/unfiltered/metrics_hs.tsv \
				   $(foreach sample,$(SAMPLES),metrics/duplex/$(sample).idx_stats.txt) \
				   $(foreach sample,$(SAMPLES),metrics/duplex/$(sample).aln_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/duplex/$(sample).insert_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/duplex/$(sample).probe-A.hs_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/duplex/$(sample).probe-B.hs_metrics.txt) \
				   metrics/duplex/metrics_idx.tsv \
				   metrics/duplex/metrics_aln.tsv \
				   metrics/duplex/metrics_insert.tsv \
				   metrics/duplex/metrics_insert_distribution.tsv \
				   metrics/duplex/metrics_hs.tsv \
				   $(foreach sample,$(SAMPLES),metrics/simplex/$(sample).idx_stats.txt) \
				   $(foreach sample,$(SAMPLES),metrics/simplex/$(sample).aln_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/simplex/$(sample).insert_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/simplex/$(sample).probe-A.hs_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/simplex/$(sample).probe-B.hs_metrics.txt) \
				   metrics/simplex/metrics_idx.tsv \
				   metrics/simplex/metrics_aln.tsv \
				   metrics/simplex/metrics_insert.tsv \
				   metrics/simplex/metrics_insert_distribution.tsv \
				   metrics/simplex/metrics_hs.tsv \
				   metrics/summary/metrics_idx.tsv \
				   metrics/summary/metrics_aln.tsv \
				   metrics/summary/metrics_insert.tsv \
				   metrics/summary/metrics_insert_distribution.tsv \
				   metrics/summary/metrics_hs.tsv \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).A.ontarget.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).B.ontarget.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).AB.offtarget.txt) \
				   metrics/summary/metrics_ts.tsv \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample)-pileup.txt) \
				   $(foreach sample,$(SAMPLES),metrics/simplex/$(sample)-pileup.txt) \
				   $(foreach sample,$(SAMPLES),metrics/duplex/$(sample)-pileup.txt)


define pileup-metric
metrics/standard/$1-pileup.txt : bam/$1-standard.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
								 cp bam/$1-standard.bam metrics/standard/$1.bam && \
								 cp bam/$1-standard.bam.bai metrics/standard/$1.bam.bai && \
								 cp bam/$1-standard.bai metrics/standard/$1.bai && \
								 cd metrics/standard && \
								 $(HOME)/share/usr/jdk1.8.0_74/bin/java -server -Xms2G -Xmx8G -cp $(WALTZ) org.mskcc.juber.waltz.Waltz PileupMetrics 20 $1.bam $(REF_FASTA) $(WALTZ_BED_FILE) && \
								 rm -rf $1.bam && \
								 rm -rf $1.bam.bai && \
								 rm -rf $1.bai && \
								 cd ..")
									 
metrics/simplex/$1-pileup.txt : bam/$1-simplex.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
								 cp bam/$1-simplex.bam metrics/simplex/$1.bam && \
								 cp bam/$1-simplex.bam.bai metrics/simplex/$1.bam.bai && \
								 cp bam/$1-simplex.bai metrics/simplex/$1.bai && \
								 cd metrics/simplex && \
								 $(HOME)/share/usr/jdk1.8.0_74/bin/java -server -Xms2G -Xmx8G -cp $(WALTZ) org.mskcc.juber.waltz.Waltz PileupMetrics 20 $1.bam $(REF_FASTA) $(WALTZ_BED_FILE) && \
								 rm -rf $1.bam && \
								 rm -rf $1.bam.bai && \
								 rm -rf $1.bai && \
								 cd ..")
								 
metrics/duplex/$1-pileup.txt : bam/$1-duplex.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
								 cp bam/$1-duplex.bam metrics/duplex/$1.bam && \
								 cp bam/$1-duplex.bam.bai metrics/duplex/$1.bam.bai && \
								 cp bam/$1-duplex.bai metrics/duplex/$1.bai && \
								 cd metrics/duplex && \
								 $(HOME)/share/usr/jdk1.8.0_74/bin/java -server -Xms2G -Xmx8G -cp $(WALTZ) org.mskcc.juber.waltz.Waltz PileupMetrics 20 $1.bam $(REF_FASTA) $(WALTZ_BED_FILE) && \
								 rm -rf $1.bam && \
								 rm -rf $1.bam.bai && \
								 rm -rf $1.bai && \
								 cd ..")
								 
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call pileup-metric,$(sample))))
				   
define coverage-metric
metrics/standard/$1.A.ontarget.txt : marianas/$1/$1.realn.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
								 samtools view -L $$(ONTARGET_FILE_A) $$(<) -b > metrics/standard/$1-ontarget-A.bam && \
								 samtools index metrics/standard/$1-ontarget-A.bam && \
								 java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) BamIndexStats \
								 I=metrics/standard/$1-ontarget-A.bam \
								 TMP_DIR=$(TMPDIR) \
								 > $$(@) && \
								 rm -rf metrics/standard/$1-ontarget-A.bam && \
								 rm -rf metrics/standard/$1-ontarget-A.bam.bai")
									 
metrics/standard/$1.B.ontarget.txt : marianas/$1/$1.realn.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
								 samtools view -L $$(ONTARGET_FILE_B) $$(<) -b > metrics/standard/$1-ontarget-B.bam && \
								 samtools index metrics/standard/$1-ontarget-B.bam && \
								 java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) BamIndexStats \
								 I=metrics/standard/$1-ontarget-B.bam \
								 TMP_DIR=$(TMPDIR) \
								 > $$(@) && \
								 rm -rf metrics/standard/$1-ontarget-B.bam && \
								 rm -rf metrics/standard/$1-ontarget-B.bam.bai")
	
metrics/standard/$1.AB.offtarget.txt : marianas/$1/$1.realn.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
								 samtools view -L $$(OFFTARGET_FILE) $$(<) -b > metrics/standard/$1-offtarget-AB.bam && \
								 samtools index metrics/standard/$1-offtarget-AB.bam && \
								 java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) BamIndexStats \
								 I=metrics/standard/$1-offtarget-AB.bam \
								 TMP_DIR=$(TMPDIR) \
								 > $$(@) && \
								 rm -rf metrics/standard/$1-offtarget-AB.bam && \
								 rm -rf metrics/standard/$1-offtarget-AB.bam.bai")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call coverage-metric,$(sample))))

define picard-metrics-standard
metrics/standard/$1.idx_stats.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) BamIndexStats \
									   I=$$(<) \
									   TMP_DIR=$(TMPDIR) \
									   > $$(@)")
									   
metrics/standard/$1.aln_metrics.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectAlignmentSummaryMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   TMP_DIR=$(TMPDIR)")

metrics/standard/$1.insert_metrics.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectInsertSizeMetrics \
									   I=$$(<) \
									   O=$$(@) \
									   H=metrics/standard/$1.insert_metrics.pdf \
									   M=0.5 \
									   TMP_DIR=$(TMPDIR)")
												
metrics/standard/$1.oxog_metrics.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectOxoGMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   TMP_DIR=$(TMPDIR)")

metrics/standard/$1.probe-A.hs_metrics.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CalculateHsMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   BAIT_INTERVALS=$(POOL_A_TARGET_FILE) \
									   TARGET_INTERVALS=$(POOL_A_TARGET_FILE) \
									   TMP_DIR=$(TMPDIR)")
												
metrics/standard/$1.probe-B.hs_metrics.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx12G -jar $$(PICARD_JAR) CalculateHsMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   BAIT_INTERVALS=$(POOL_B_TARGET_FILE) \
									   TARGET_INTERVALS=$(POOL_B_TARGET_FILE) \
									   TMP_DIR=$(TMPDIR)")
									   
metrics/standard/$1.probe-A.hs_metrics-nodedup.txt : marianas/$1/$1.realn.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CalculateHsMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   BAIT_INTERVALS=$(POOL_A_TARGET_FILE) \
									   TARGET_INTERVALS=$(POOL_A_TARGET_FILE) \
									   TMP_DIR=$(TMPDIR)")
												
metrics/standard/$1.probe-B.hs_metrics-nodedup.txt : marianas/$1/$1.realn.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx12G -jar $$(PICARD_JAR) CalculateHsMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   BAIT_INTERVALS=$(POOL_B_TARGET_FILE) \
									   TARGET_INTERVALS=$(POOL_B_TARGET_FILE) \
									   TMP_DIR=$(TMPDIR)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call picard-metrics-standard,$(sample))))
		
define picard-metrics-unfiltered
metrics/unfiltered/$1.idx_stats.txt : bam/$1-unfiltered.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) BamIndexStats \
									   I=$$(<) \
									   TMP_DIR=$(TMPDIR) \
									   > $$(@)")
									   
metrics/unfiltered/$1.aln_metrics.txt : bam/$1-unfiltered.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectAlignmentSummaryMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   TMP_DIR=$(TMPDIR)")

metrics/unfiltered/$1.insert_metrics.txt : bam/$1-unfiltered.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectInsertSizeMetrics \
									   I=$$(<) \
									   O=$$(@) \
									   H=metrics/unfiltered/$1.insert_metrics.pdf \
									   M=0.5 \
									   TMP_DIR=$(TMPDIR)")
												
metrics/unfiltered/$1.probe-A.hs_metrics.txt : bam/$1-unfiltered.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CalculateHsMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   BAIT_INTERVALS=$(POOL_A_TARGET_FILE) \
									   TARGET_INTERVALS=$(POOL_A_TARGET_FILE) \
									   TMP_DIR=$(TMPDIR)")
												
metrics/unfiltered/$1.probe-B.hs_metrics.txt : bam/$1-unfiltered.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx12G -jar $$(PICARD_JAR) CalculateHsMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   BAIT_INTERVALS=$(POOL_B_TARGET_FILE) \
									   TARGET_INTERVALS=$(POOL_B_TARGET_FILE) \
									   TMP_DIR=$(TMPDIR)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call picard-metrics-unfiltered,$(sample))))		

define picard-metrics-duplex
metrics/duplex/$1.idx_stats.txt : bam/$1-duplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) BamIndexStats \
									   I=$$(<) \
									   TMP_DIR=$(TMPDIR) \
									   > $$(@)")
									   
metrics/duplex/$1.aln_metrics.txt : bam/$1-duplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectAlignmentSummaryMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   TMP_DIR=$(TMPDIR)")

metrics/duplex/$1.insert_metrics.txt : bam/$1-duplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectInsertSizeMetrics \
									   I=$$(<) \
									   O=$$(@) \
									   H=metrics/duplex/$1.insert_metrics.pdf \
									   M=0.5 \
									   TMP_DIR=$(TMPDIR)")
												
metrics/duplex/$1.probe-A.hs_metrics.txt : bam/$1-duplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CalculateHsMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   BAIT_INTERVALS=$(POOL_A_TARGET_FILE) \
									   TARGET_INTERVALS=$(POOL_A_TARGET_FILE) \
									   TMP_DIR=$(TMPDIR)")
												
metrics/duplex/$1.probe-B.hs_metrics.txt : bam/$1-duplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx12G -jar $$(PICARD_JAR) CalculateHsMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   BAIT_INTERVALS=$(POOL_B_TARGET_FILE) \
									   TARGET_INTERVALS=$(POOL_B_TARGET_FILE) \
									   TMP_DIR=$(TMPDIR)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call picard-metrics-duplex,$(sample))))
		
define picard-metrics-simplex
metrics/simplex/$1.idx_stats.txt : bam/$1-simplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) BamIndexStats \
									   I=$$(<) \
									   TMP_DIR=$(TMPDIR) \
									   > $$(@)")
									   
metrics/simplex/$1.aln_metrics.txt : bam/$1-simplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectAlignmentSummaryMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   TMP_DIR=$(TMPDIR)")

metrics/simplex/$1.insert_metrics.txt : bam/$1-simplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectInsertSizeMetrics \
									   I=$$(<) \
									   O=$$(@) \
									   H=metrics/simplex/$1.insert_metrics.pdf \
									   M=0.5 \
									   TMP_DIR=$(TMPDIR)")
												
metrics/simplex/$1.probe-A.hs_metrics.txt : bam/$1-simplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CalculateHsMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   BAIT_INTERVALS=$(POOL_A_TARGET_FILE) \
									   TARGET_INTERVALS=$(POOL_A_TARGET_FILE) \
									   TMP_DIR=$(TMPDIR)")
												
metrics/simplex/$1.probe-B.hs_metrics.txt : bam/$1-simplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx12G -jar $$(PICARD_JAR) CalculateHsMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   BAIT_INTERVALS=$(POOL_B_TARGET_FILE) \
									   TARGET_INTERVALS=$(POOL_B_TARGET_FILE) \
									   TMP_DIR=$(TMPDIR)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call picard-metrics-simplex,$(sample))))
		
metrics/standard/metrics_idx.tsv : $(wildcard metrics/standard/$(SAMPLES).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 1 --sample_names '$(SAMPLES)'")
		
metrics/standard/metrics_aln.tsv : $(wildcard metrics/standard/$(SAMPLES).aln_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 2 --sample_names '$(SAMPLES)'")
	
metrics/standard/metrics_insert.tsv : $(wildcard metrics/standard/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 3 --sample_names '$(SAMPLES)'")
	
metrics/standard/metrics_insert_distribution.tsv : $(wildcard metrics/standard/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 4 --sample_names '$(SAMPLES)'")
	
metrics/standard/metrics_oxog.tsv : $(wildcard metrics/standard/$(SAMPLES).oxog_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 5 --sample_names '$(SAMPLES)'")

metrics/standard/metrics_hs.tsv : $(wildcard metrics/standard/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/standard/$(SAMPLES).probe-B.hs_metrics.txt) $(wildcard metrics/standard/$(SAMPLES).probe-A.hs_metrics-nodedup.txt) $(wildcard metrics/standard/$(SAMPLES).probe-B.hs_metrics-nodedup.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 6 --sample_names '$(SAMPLES)'")
	
metrics/unfiltered/metrics_idx.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 7 --sample_names '$(SAMPLES)'")
		
metrics/unfiltered/metrics_aln.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).aln_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 8 --sample_names '$(SAMPLES)'")
	
metrics/unfiltered/metrics_insert.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 9 --sample_names '$(SAMPLES)'")
	
metrics/unfiltered/metrics_insert_distribution.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 10 --sample_names '$(SAMPLES)'")
	
metrics/unfiltered/metrics_hs.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/unfiltered/$(SAMPLES).probe-B.hs_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 11 --sample_names '$(SAMPLES)'")
	
metrics/duplex/metrics_idx.tsv : $(wildcard metrics/duplex/$(SAMPLES).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 12 --sample_names '$(SAMPLES)'")
		
metrics/duplex/metrics_aln.tsv : $(wildcard metrics/duplex/$(SAMPLES).aln_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 13 --sample_names '$(SAMPLES)'")
	
metrics/duplex/metrics_insert.tsv : $(wildcard metrics/duplex/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 14 --sample_names '$(SAMPLES)'")
	
metrics/duplex/metrics_insert_distribution.tsv : $(wildcard metrics/duplex/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 15 --sample_names '$(SAMPLES)'")
	
metrics/duplex/metrics_hs.tsv : $(wildcard metrics/duplex/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/duplex/$(SAMPLES).probe-B.hs_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 16 --sample_names '$(SAMPLES)'")
	
metrics/simplex/metrics_idx.tsv : $(wildcard metrics/simplex/$(SAMPLES).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 17 --sample_names '$(SAMPLES)'")
		
metrics/simplex/metrics_aln.tsv : $(wildcard metrics/simplex/$(SAMPLES).aln_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 18 --sample_names '$(SAMPLES)'")
	
metrics/simplex/metrics_insert.tsv : $(wildcard metrics/simplex/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 19 --sample_names '$(SAMPLES)'")
	
metrics/simplex/metrics_insert_distribution.tsv : $(wildcard metrics/simplex/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 20 --sample_names '$(SAMPLES)'")
	
metrics/simplex/metrics_hs.tsv : $(wildcard metrics/simplex/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/simplex/$(SAMPLES).probe-B.hs_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 21 --sample_names '$(SAMPLES)'")
	
metrics/summary/metrics_idx.tsv : metrics/standard/metrics_idx.tsv metrics/unfiltered/metrics_idx.tsv metrics/duplex/metrics_idx.tsv metrics/simplex/metrics_idx.tsv
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 22")
		
metrics/summary/metrics_aln.tsv : metrics/standard/metrics_aln.tsv metrics/unfiltered/metrics_aln.tsv metrics/duplex/metrics_aln.tsv metrics/simplex/metrics_aln.tsv
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 23")
	
metrics/summary/metrics_insert.tsv : metrics/standard/metrics_insert.tsv metrics/unfiltered/metrics_insert.tsv metrics/duplex/metrics_insert.tsv metrics/simplex/metrics_insert.tsv
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 24")
	
metrics/summary/metrics_insert_distribution.tsv : metrics/standard/metrics_insert_distribution.tsv metrics/unfiltered/metrics_insert_distribution.tsv metrics/duplex/metrics_insert_distribution.tsv metrics/simplex/metrics_insert_distribution.tsv
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 25")
	
metrics/summary/metrics_hs.tsv : metrics/standard/metrics_hs.tsv metrics/unfiltered/metrics_hs.tsv metrics/duplex/metrics_hs.tsv metrics/simplex/metrics_hs.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 26")

metrics/summary/metrics_ts.tsv : $(wildcard metrics/standard/$(SAMPLES).A.ontarget.txt) $(wildcard metrics/standard/$(SAMPLES).B.ontarget.txt) $(wildcard metrics/standard/$(SAMPLES).AB.offtarget.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 27 --sample_names '$(SAMPLES)'")
	
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
