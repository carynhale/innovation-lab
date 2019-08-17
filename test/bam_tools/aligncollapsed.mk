include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/align_collapsed.$(NOW)
PHONY += marianas

BWA_ALN_OPTS ?= -M
BWAMEM_REF_FASTA ?= $(REF_FASTA)
BWAMEM_THREADS = 12
BWAMEM_MEM_PER_THREAD = 2G

SAMTOOLS_THREADS = 8
SAMTOOLS_MEM_THREAD = 2G

GATK_THREADS = 8
GATK_MEM_THREAD = 2G

JAVA = /home/${USER}/share/usr/jdk1.8.0_74/bin/java
BED_FILE=/home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-AB.waltz.bed

align_collapsed : $(foreach sample,$(SAMPLES),marianas/$(sample)/timestamp)

define fastq-to-bam
marianas/$1/$1.collapsed.bwamem.bam : marianas/$1/collapsed_R1_.fastq marianas/$1/collapsed_R1_.fastq
	$$(call RUN,-c -n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD) -w 1440,"set -o pipefail && \
																		           $(BWA) mem -t $(BWAMEM_THREADS) $(BWA_ALN_OPTS) -R \"@RG\tID:$1\tLB:$1\tPL:${SEQ_PLATFORM}\tSM:$1\" $(BWAMEM_REF_FASTA) $$(^) | $(SAMTOOLS) view -bhS - > $$(@)")
																		           
marianas/$1/$1.collapsed.sorted.bam : marianas/$1/$1.collapsed.bwamem.bam
	$$(call RUN,-c -n $(SAMTOOLS_THREADS) -s 1G -m $(SAMTOOLS_MEM_THREAD) -w 1440,"set -o pipefail && \
								  									   		   	   samtools sort -@ $(SAMTOOLS_THREADS) -m $(SAMTOOLS_MEM_THREAD) $$(^) -o $$(@) -T $(TMPDIR) && \
								  									   		   	   samtools index $$(@) && \
								  									   		   	   cp marianas/$1/$1.collapsed.sorted.bam.bai marianas/$1/$1.collapsed.sorted.bai")
								  									   		   
marianas/$1/$1.collapsed.fixed.bam : marianas/$1/$1.collapsed.sorted.bam
	$$(call RUN,-c -n 1 -s 12G -m 18G -w 1440,"set -o pipefail && \
									  		   /home/${USER}/share/usr/jdk1.8.0_74/bin/java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx16G -jar /home/${USER}/share/usr/picard/bin/picard.jar \
									  		   FixMateInformation \
									  		   I=$$(<) \
									  		   O=$$(@) \
									  		   SORT_ORDER=coordinate \
									  		   TMP_DIR=$(TMPDIR) \
									  		   COMPRESSION_LEVEL=0 \
									  		   CREATE_INDEX=true \
									  		   VALIDATION_STRINGENCY=LENIENT")

marianas/$1/timestamp : marianas/$1/$1.collapsed.fixed.bam
	$$(call RUN,-c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  cd marianas/$1 && \
									  $(JAVA) -server -Xms2G -Xmx8G -cp $(MARIANAS) org.mskcc.marianas.umi.duplex.postprocessing.SeparateBams $1.collapsed.fixed.bam && \
									  echo 'Done!\n' > timestamp && \
									  cd ../..")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fastq-to-bam,$(sample))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
