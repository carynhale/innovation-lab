include modules/Makefile.inc

LOGDIR ?= log/mspyclone.$(NOW)
PHONY += pyclone sufam

MIN_DEPTH ?= 50

pyclone : $(foreach set,$(SAMPLE_SETS),sufam/$(set).tsv) \
		  $(foreach set,$(SAMPLE_SETS),pyclone/$(set)/config.yaml) \
		  $(foreach set,$(SAMPLE_SETS),pyclone/$(set)/trace/alpha.tsv.bz2) \
		  $(foreach set,$(SAMPLE_SETS),pyclone/$(set)/report/pyclone.tsv) \
		  $(foreach set,$(SAMPLE_SETS),pyclone/$(set)/report/pyclone.pdf)
		  
define multi-sample-pyclone
sufam/%.txt : summary/tsv/mutation_summary.tsv
	$$(call RUN,-s 8G -m 16G,"$(RSCRIPT) $(SCRIPTS_DIR)/variant_callers/combinesamples.R --sample_set $$*")

sufam/%.tsv : sufam/%.txt
	$$(call RUN,-s 8G -m 16G,"$(RSCRIPT) $(SCRIPTS_DIR)/variant_callers/updatesamples.R --sample_set $$*")
	
pyclone/%/config.yaml : sufam/%.tsv
	$$(call RUN, -s 4G -m 6G,"mkdir -p pyclone/$$(*) && \
							  $(RSCRIPT) $(SCRIPTS_DIR)/clonality/config2mspyclone.R --sample_set $$(*) --normal_samples $(NORMAL_SAMPLES) && \
							  $(RSCRIPT) $(SCRIPTS_DIR)/clonality/tsvforpyclone.R --sample_set $$(*) --normal_samples $(NORMAL_SAMPLES) --min_depth $(MIN_DEPTH)")

pyclone/%/trace/alpha.tsv.bz2 : pyclone/%/config.yaml
	$$(call RUN,-s 4G -m 6G -w 7200,"source /home/$(USER)/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/$(USER)/share/usr/anaconda-envs/PyClone-0.13.1 && \
							 		 PyClone run_analysis --config_file pyclone/$$*/config.yaml --seed 0")

pyclone/%/report/pyclone.tsv : pyclone/%/trace/alpha.tsv.bz2
	$$(call RUN,-s 4G -m 6G -w 7200,"make -p pyclone/$$*/report && \
									 source /home/$(USER)/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/$(USER)/share/usr/anaconda-envs/PyClone-0.13.1 && \
							 		 PyClone build_table --config_file pyclone/$$*/config.yaml --out_file pyclone/$$*/report/pyclone.tsv --max_cluster 10 --table_type old_style --burnin 5000")
							 		 
pyclone/%/report/pyclone.pdf : pyclone/%/report/pyclone.tsv
	$$(call RUN,-s 4G -m 6G -w 7200,"$(RSCRIPT) $(SCRIPTS_DIR)/clonality/plot2mspyclone.R --sample_set $$(*) --normal_samples $(NORMAL_SAMPLES) --min_depth $(MIN_DEPTH)")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call multi-sample-pyclone,$(set))))
		

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
