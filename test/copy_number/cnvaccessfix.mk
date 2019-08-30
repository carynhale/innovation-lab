include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_fix.$(NOW)
PHONY += cnvaccess cnvaccess/cnr

cnvaccess_fix : $(foreach sample,$(SAMPLES),cnvaccess/cnr/$(sample).pool-A.cnr) \
				$(foreach sample,$(SAMPLES),cnvaccess/cnr/$(sample).pool-B.cnr) \
				$(foreach sample,$(SAMPLES),cnvaccess/cnr/$(sample).no-pool.cnr)
				
ACCESS_REF_FILE_A = ~/share/reference/cnvkit_reference/MSK-ACCESS-v1_0-probe-A.cnr
ACCESS_REF_FILE_B = ~/share/reference/cnvkit_reference/MSK-ACCESS-v1_0-probe-B.cnr
ACCESS_REF_FILE_OFF = ~/share/reference/cnvkit_reference/MSK-ACCESS-v1_0-noprobe.cnr

define cnvaccess-cnr
cnvaccess/cnr/%.pool-A.cnr : cnvaccess/cnn/%.pool-A.targetcoverage.cnn cnvaccess/cnn/%.pool-A.antitargetcoverage.cnn
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) $(ACCESS_REF_FILE_A) -o cnvaccess/cnr/$$(*).pool-A.cnr")

cnvaccess/cnr/%.pool-B.cnr : cnvaccess/cnn/%.pool-B.targetcoverage.cnn cnvaccess/cnn/%.pool-B.antitargetcoverage.cnn
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) $(ACCESS_REF_FILE_B) -o cnvaccess/cnr/$$(*).pool-B.cnr")
	
cnvaccess/cnr/%.no-pool.cnr : cnvaccess/cnn/%.no-pool.targetcoverage.cnn cnvaccess/cnn/%.no-pool.antitargetcoverage.cnn
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) $(ACCESS_REF_FILE_OFF) -o cnvaccess/cnr/$$(*).no-pool.cnr")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvaccess-cnr,$(sample))))
				
.PHONY: $(PHONY)
