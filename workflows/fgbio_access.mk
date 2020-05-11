include innovation-lab/Makefile.inc
include innovation-lab/config/fgbio.inc
include innovation-lab/config/gatk.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/fgbio_access.$(NOW)

fgbio_access : $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_R1.fastq.gz)

define copy-fastq
fgbio/$1/$1_R1.fastq.gz : $3
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
									 mkdir -p fgbio/$1 && \
								     $(RSCRIPT) $(SCRIPTS_DIR)/fastq_tools/copy_fastq.R --sample_name $1 --directory_name 'fgbio' --fastq_files '$$^'")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call copy-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))

define fastq-2bam
fgbio/$1/ : fgbio/$1/$1_R1.fastq.gz
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $$(call FGBIO_CMD,2G,8G) \
									  FastqToBam \
									  --input fgbio/$1/$1_R1.fastq.gz fgbio/$1/$1_R2.fastq.gz \
									  --read-structures 3M2S+T 3M2S+T \
									  --output ${OUTPUT_FILE} \
									  --sample ${SAMPLE_NAME} \
									  --library ${SAMPLE_NAME} \
									  --platform ${sequencing_platform} \
									  --platform-unit ${sequencing_platform_unit} \
									  --sequencing-center ${sequencing_center} \
									  --read-group-id ${READ_GROUP}")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call clip-umi,$(sample))))
	
COMMAND="$COMMAND \"$JAVA -Djava.io.tmpdir=/scratch/ -Xmx${MEM}g -jar
${FGBIO} FastqToBam
--input ${R1_FASTQ} ${R2_FASTQ}

\""


..DUMMY := $(shell mkdir -p version)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: fgbio_access
