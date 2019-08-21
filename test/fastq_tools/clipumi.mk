include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/clip_umi.$(NOW)
PHONY += marianas

clip_umi : $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)_R1.fastq.gz) \
		   $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)_R2.fastq.gz) \
		   $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)_R1_umi-clipped.fastq.gz) \
		   $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)_R2_umi-clipped.fastq.gz)

JAVA = $(HOME)/share/usr/jdk1.8.0_74/bin/java
MARIANAS_UMI_LENGTH ?= 3
MARIANAS_MIN_MAPQ ?= 1
MARIANAS_MIN_BAQ ?= 20
MARIANAS_MISMATCH ?= 0
MARIANAS_WOBBLE ?= 1
MARIANAS_MIN_CONSENSUS ?= 90
WALTZ_MIN_MAPQ ?= 20

define copy-fastq
marianas/$1/$1_R1.fastq.gz marianas/$1/$1_R2.fastq.gz : $3
	$$(call RUN,-c -n 1 -s 2G -m 4G,"mkdir -p marianas/$1 && \
								     $(RSCRIPT) modules/test/fastq_tools/copyfastq.R --sample_name $1 --fastq_files '$$^'")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call copy-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))
 

define clip-umi
marianas/$1/$1_R1_umi-clipped.fastq.gz marianas/$1/$1_R2_umi-clipped.fastq.gz : marianas/$1/$1_R1.fastq.gz marianas/$1/$1_R2.fastq.gz
	$$(call RUN,-c -n 1 -s 8G -m 16G,"cd marianas/$1/ && \
									  $(JAVA) -Djava.io.tmpdir=$(TMPDIR) -server -Xms2G -Xmx8G -cp $(MARIANAS) \
									  org.mskcc.marianas.umi.duplex.fastqprocessing.ProcessLoopUMIFastq \
									  $1_R1.fastq.gz $1_R2.fastq.gz \
									  $(MARIANAS_UMI_LENGTH) && \
									  cd ../..")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call clip-umi,$(sample))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
