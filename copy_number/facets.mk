include innovation-lab/Makefile.inc

LOGDIR ?= log/facets.$(NOW)

facets : facets/targets/targets.vcf

facets/targets/targets.vcf : $(TARGETS_FILE)
	$(INIT) $(BEDTOOLS) intersect -header -u -a $(DBSNP) -b $< > $@
		 