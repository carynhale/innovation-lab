include innovation-lab/Makefile.inc

LOGDIR ?= log/link_bam.$(NOW)

link_bam : $(foreach sample,$(SAMPLES),bam/$(sample).bam)

define link-bam
bam/%.bam : /ifs/dmpshare/share/irb12_245/%.bam
	$$(call RUN,-c -n 1 -s 1G -m 2G,"ln -s /ifs/dmpshare/share/irb12_245/$$(*).bam bam/$$(*).bam && \
									 ln -s /ifs/dmpshare/share/irb12_245/$$(*).bai bam/$$(*).bai")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call link-bam,$(sample))))


..DUMMY := $(shell mkdir -p version; \
			 ln --version > version/link_bam.txt)
.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: link_bam
