include innovation-lab/Makefile.inc




..DUMMY := $(shell mkdir -p version; \
			 $(BWA) &> version/tmp.txt; \
			 head -3 version/tmp.txt | tail -2 > version/msk_access.txt; \
			 rm version/tmp.txt; \
			 $(SAMTOOLS) --version >> version/msk_access.txt; \
			 echo "gatk3" >> version/msk_access.txt; \
			 $(GATK) --version >> version/msk_access.txt; \
			 echo "picard" >> version/msk_access.txt; \
			 R --version >> version.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: msk_access
