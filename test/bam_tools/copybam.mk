include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/copy_bam.$(NOW)
PHONY += bam

JAVA = /home/${USER}/share/usr/jdk1.8.0_74/bin/java

copy_bam : $(foreach sample,$(SAMPLES),bam/$(sample)-standard.bam) \
		   $(foreach sample,$(SAMPLES),bam/$(sample)-unfiltered.bam)
		   $(foreach sample,$(SAMPLES),bam/$(sample)-simplex.bam) \
		   $(foreach sample,$(SAMPLES),bam/$(sample)-duplex.bam)


define copy-to-bam
bam/$1-standard.bam : marianas/$1/$1.standard.bam
	$$(call RUN, -c -n 1 -s 12G -m 18G -w 1440,"set -o pipefail && \
												java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx16G -jar $$(PICARD_JAR) AddOrReplaceReadGroups \
												I=$$(<) \
												O=$$(@) \
												RGID=$1 \
												RGLB=$1 \
												RGPL=illumina \
												RGPU=NA \
												RGSM=$1 \
												TMP_DIR=$(TMPDIR) && \
												samtools index $$(@) && \
												cp bam/$1-standard.bam.bai bam/$1-standard.bai")

bam/$1-unfiltered.bam : marianas/$1/$1.collapsed.bam
	$$(call RUN, -c -s 2G -m 4G,"set -o pipefail && \
								 cp $$(<) $$(@) && \
								 cp marianas/$1/$1.collapsed.bam.bai bam/$1-unfiltered.bam.bai && \
								 cp marianas/$1/$1.collapsed.bai bam/$1-unfiltered.bai")
												
bam/$1-simplex.bam : marianas/$1/timestamp
	$$(call RUN, -c -n 1 -s 12G -m 18G -w 1440,"set -o pipefail && \
												java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx16G -jar $$(PICARD_JAR) AddOrReplaceReadGroups \
												I=marianas/$1/$1.collapsed-simplex.bam \
												O=$$(@) \
												RGID=$1 \
												RGLB=$1 \
												RGPL=illumina \
												RGPU=NA \
												RGSM=$1 \
												TMP_DIR=$(TMPDIR) && \
												samtools index $$(@) && \
												cp bam/$1-simplex.bam.bai bam/$1-simplex.bai")
												
bam/$1-duplex.bam : marianas/$1/timestamp
	$$(call RUN, -c -n 1 -s 12G -m 18G -w 1440,"set -o pipefail && \
												java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx16G -jar $$(PICARD_JAR) AddOrReplaceReadGroups \
												I=marianas/$1/$1.collapsed-duplex.bam \
												O=$$(@) \
												RGID=$1 \
												RGLB=$1 \
												RGPL=illumina \
												RGPU=NA \
												RGSM=$1 \
												TMP_DIR=$(TMPDIR) && \
												samtools index $$(@) && \
												cp bam/$1-duplex.bam.bai bam/$1-duplex.bai")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call copy-to-bam,$(sample))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
