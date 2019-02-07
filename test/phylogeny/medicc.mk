include modules/Makefile.inc

INDEX = $(seq -f "%03g" 1 100)

LOGDIR ?= log/medicc.$(NOW)
PHONY += medicc medicc/mad medicc/mpcf medicc/medicc medicc/boot

medicc : $(foreach set,$(SAMPLE_SETS),medicc/mad/$(set).RData) $(foreach set,$(SAMPLE_SETS),medicc/mpcf/$(set).RData) $(foreach set,$(SAMPLE_SETS),medicc/medicc/$(set)/desc.txt) $(foreach set,$(SAMPLE_SETS),medicc/medicc/$(set)/tree_final.new) $(foreach set,$(SAMPLE_SETS),medicc/boot/$(set)/$(INDEX)) $(foreach set,$(SAMPLE_SETS),medicc/boot/$(set)/$INDEX/tree_final.new)


define combine-samples
medicc/mad/%.RData : $(wildcard $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).Rdata))
	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"$(RSCRIPT) modules/test/phylogeny/combinesamples.R --sample_set $$* --normal_samples '$(NORMAL_SAMPLES)'")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call combine-samples,$(set))))

define ascat-mpcf
medicc/mpcf/%.RData : medicc/mad/%.RData
	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"if [ ! -d medicc/ascat ]; then mkdir medicc/ascat; fi && \
												 $(RSCRIPT) modules/test/phylogeny/segmentsamples.R --sample_set $$* --normal_samples '$(NORMAL_SAMPLES)' --gamma '$${mpcf_gamma}' --nlog2 '$${mpcf_nlog2}' --nbaf '$${mpcf_nbaf}'")
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call ascat-mpcf,$(set))))

define init-medicc
medicc/medicc/%/desc.txt : medicc/mpcf/%.RData
	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"$(RSCRIPT) modules/test/phylogeny/initmedicc.R --sample_set $$*")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call init-medicc,$(set))))

define run-medicc
medicc/medicc/%/tree_final.new : medicc/medicc/%/desc.txt
	$$(call RUN,-c -s 8G -m 12G -v $(MEDICC_ENV),"source $(MEDICC_VAR) && \
												  $(MEDICC_BIN)/medicc.py medicc/medicc/$$*/desc.txt medicc/medicc/$$* -v")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call run-medicc,$(set))))

define boot-medicc
$(wildcard medicc/boot/%/$(INDEX)) : medicc/mpcf/%.RData
	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"$(RSCRIPT) modules/test/phylogeny/bootstrapmedicc.R --sample_set $$*")

$(wildcard medicc/boot/%/$(INDEX)/tree_final.new) : $(wildcard medicc/boot/%/$(INDEX))
	$$(call RUN,-c -s 2G -m 4G -n 12 -v $(MEDICC_ENV),"source $(MEDICC_VAR) && \
												  	   seq -f '%03g' 1 100 | parallel -j 12 'if [ ! -f medicc/boot/$$*/{}/tree_final.new ]; then $(MEDICC_BIN)/medicc.py medicc/boot/$$*/{}/desc.txt medicc/boot/$$*/{}/ -v; fi'")
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call boot-medicc,$(set))))
		
		
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)





  