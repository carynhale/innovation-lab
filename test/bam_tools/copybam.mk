include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/copy_bam.$(NOW)
PHONY += bam

copy_bam : $(foreach sample,$(SAMPLES),bam/$(sample)-standard.bam) \
		   $(foreach sample,$(SAMPLES),bam/$(sample)-unfiltered.bam)
		   $(foreach sample,$(SAMPLES),bam/$(sample)-simplex.bam) \
		   $(foreach sample,$(SAMPLES),bam/$(sample)-duplex.bam)


define copy-to-bam
bam/$1-standard.bam : marianas/$1/$1.standard.bam
	$$(call RUN, -c -s 2G -m 4G ,"set -o pipefail && \
								 cp $$(<) $$(@) && \
								 cp marianas/$1/$1.standard.bam.bai bam/$1-standard.bam.bai && \
								 cp marianas/$1/$1.standard.bai bam/$1-standard.bai")

bam/$1-unfiltered.bam : marianas/$1/$1.collapsed.bam
	$$(call RUN, -c -s 2G -m 4G,"set -o pipefail && \
								 cp $$(<) $$(@) && \
								 cp marianas/$1/$1.collapsed.bam.bai bam/$1-unfiltered.bam.bai && \
								 cp marianas/$1/$1.collapsed.bai bam/$1-unfiltered.bai")
												
bam/$1-simplex.bam : marianas/$1/timestamp
	$$(call RUN, -c -s 2G -m 4G,"set -o pipefail && \
								 cp marianas/$1/$1.collapsed-simplex $$(@) && \
								 cp marianas/$1/$1.collapsed-simplex.bai bam/$1-simplex.bam.bai && \
								 cp marianas/$1/$1.collapsed-simplex.bai bam/$1-simplex.bai")
												
bam/$1-duplex.bam : marianas/$1/timestamp
	$$(call RUN, -c -s 2G -m 4G,"set -o pipefail && \
								 cp marianas/$1/$1.collapsed-duplex $$(@) && \
								 cp marianas/$1/$1.collapsed-duplex.bai bam/$1-duplex.bam.bai && \
								 cp marianas/$1/$1.collapsed-duplex.bai bam/$1-duplex.bai")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call copy-to-bam,$(sample))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
