include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/bam_interval_metrics.$(NOW)
PHONY += metrics metrics/pileup metrics/summary

#POOL_A_TARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.list
#POOL_B_TARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.list
#ONTARGET_FILE_A ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.bed
#ONTARGET_FILE_B ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.bed
#OFFTARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-AB.offtarget.bed

interval_metrics : $(foreach sample,$(SAMPLES),metrics/standard/$(sample)-pileup.txt)
#				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).idx_stats.txt) \
#				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).aln_metrics.txt) \
#				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).insert_metrics.txt) \
#				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).oxog_metrics.txt) \
#				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).probe-A.hs_metrics.txt) \
#				   metrics/standard/metrics_idx.tsv \
#				   metrics/standard/metrics_aln.tsv \
#				   metrics/standard/metrics_insert.tsv \
#				   metrics/standard/metrics_insert_distribution.tsv \
#				   metrics/standard/metrics_oxog.tsv \
#				   metrics/standard/metrics_hs.tsv \
#				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).A.ontarget.txt) \
#				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).B.ontarget.txt) \
#				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).AB.offtarget.txt) \
#				   metrics/summary/metrics_ts.tsv
				   


define pileup-metric
metrics/pielup/$1.txt : bam/$1.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
								 cp bam/$1.bam metrics/pielup/$1.bam && \
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
		
		
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)