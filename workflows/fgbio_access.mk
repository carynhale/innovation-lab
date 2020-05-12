include innovation-lab/Makefile.inc
include innovation-lab/config/fgbio.inc
include innovation-lab/config/gatk.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/fgbio_access.$(NOW)

fgbio_access : $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_R1.fastq.gz) \
			   $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_fq.bam) \
			   $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl.fastq.gz)

define copy-fastq
fgbio/$1/$1_R1.fastq.gz : $3
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
									 mkdir -p fgbio/$1 && \
								     $(RSCRIPT) $(SCRIPTS_DIR)/fastq_tools/copy_fastq.R \
								     --sample_name $1 \
								     --directory_name fgbio \
								     --fastq_files '$$^'")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call copy-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))

define fastq-2-bam
fgbio/$1/$1_fq.bam : fgbio/$1/$1_R1.fastq.gz
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $$(call FGBIO_CMD,2G,8G) \
									  FastqToBam \
									  --input fgbio/$1/$1_R1.fastq.gz fgbio/$1/$1_R2.fastq.gz \
									  --read-structures 3M2S+T 3M2S+T \
									  --output $$(@) \
									  --sample $1 \
									  --library $1 \
									  --platform illumina \
									  --platform-unit NA")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call fastq-2-bam,$(sample))))
	
define mark-adapters
fgbio/$1/$1_cl.fastq.gz : fgbio/$1/$1_fq.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $$(MARK_ADAPTERS) \
									  INPUT=$$(<) \
									  OUTPUT=/dev/stdout \
									  METRICS=fgbio/$1/$1_adapter-metrics.txt | \
									  $$(SAM_TO_FASTQ) \
									  INPUT=/dev/stdin \
									  FASTQ=$$(@) \
									  INTERLEAVE=true \
									  CLIPPING_ATTRIBUTE=XT \
									  CLIPPING_ACTION=X \
									  CLIPPING_MIN_LENGTH=25")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call mark-adapters,$(sample))))

	
..DUMMY := $(shell mkdir -p version; \
			 $(JAVA8) -jar $(FGBIO) --help &> version/fgbio_access.txt; \
			 echo "picard" >> version/fgbio_access.txt; \
			 $(PICARD) MarkIlluminaAdapters --version &>> version/fgbio_access.txt; \
			 $(PICARD) SamToFastq --version &>> version/fgbio_access.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: fgbio_access
