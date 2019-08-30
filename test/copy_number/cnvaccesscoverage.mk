include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_coverage.$(NOW)
PHONY += cnvaccess cnvaccess/cnn

cnvaccess_coverage : $(foreach sample,$(SAMPLES),cnvaccess/cnn/$(sample).pool-A.targetcoverage.cnn cnvaccess/cnn/$(sample).pool-B.targetcoverage.cnn cnvaccess/cnn/$(sample).no-pool.antitargetcoverage.cnn) \
					 $(foreach sample,$(SAMPLES),cnvaccess/cnn/$(sample).pool-A.antitargetcoverage.cnn cnvaccess/cnn/$(sample).pool-B.antitargetcoverage.cnn cnvaccess/cnn/$(sample).no-pool.targetcoverage.cnn)

ONTARGET_FILE_A ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.bed
ONTARGET_FILE_B ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.bed
OFFTARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-AB.offtarget.bed

define cnvaccess-cnn
cnvaccess/cnn/%.pool-A.targetcoverage.cnn : bam/%-standard.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE_A) -o cnvaccess/cnn/$$(*).pool-A.targetcoverage.cnn")

cnvaccess/cnn/%.pool-A.antitargetcoverage.cnn : bam/%-standard.bam
	$$(call RUN,-n 1 -s 1G -m 2G,"touch cnvaccess/cnn/$$(*).pool-A.antitargetcoverage.cnn")
	
cnvaccess/cnn/%.pool-B.targetcoverage.cnn : bam/%-standard.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE_B) -o cnvaccess/cnn/$$(*).pool-B.targetcoverage.cnn")

cnvaccess/cnn/%.pool-B.antitargetcoverage.cnn : bam/%-standard.bam
	$$(call RUN,-n 1 -s 1G -m 2G,"touch cnvaccess/cnn/$$(*).pool-B.antitargetcoverage.cnn")

cnvaccess/cnn/%.no-pool.antitargetcoverage.cnn : bam/%-standard.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(OFFTARGET_FILE) -o cnvaccess/cnn/$$(*).no-pool.antitargetcoverage.cnn")

cnvaccess/cnn/%.no-pool.targetcoverage.cnn : bam/%-standard.bam
	$$(call RUN,-n 1 -s 1G -m 2G,"touch cnvaccess/cnn/$$(*).no-pool.targetcoverage.cnn")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvaccess-cnn,$(sample))))
		
.PHONY: $(PHONY)
