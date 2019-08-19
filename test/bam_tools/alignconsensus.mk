include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/align_consensus.$(NOW)
PHONY += fgbio bam

align_consensus : $(foreach sample,$(SAMPLES),bam/$(sample).bam)

JAVA = /home/${USER}/share/usr/jdk1.8.0_74/bin/java
PICARD = /home/${USER}/share/usr/picard/bin/picard.jar
POOL_A_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.list
POOL_B_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.list

define bam-to-bam
fgbio/%.sorted.bam : fgbio/%.filtered.bam
	$$(call RUN,-c -s 12G -m 16G,"set -o pipefail && \
									  $(JAVA) -Xmx12G -jar $(PICARD) SortSam \
									  I=fgbio/$$(*).filtered.bam \
									  O=fgbio/$$(*).sorted.bam \
									  SORT_ORDER=queryname \
									  TMP_DIR=$(TMPDIR)")

fgbio/%.resorted.bam : fgbio/%.sorted.bam
	$$(call RUN,-c -n 12 -s 2G -m 4G,"set -o pipefail && \
									  $(JAVA) -Xmx16G -jar $(PICARD) SamToFastq \
									  I=fgbio/$$(*).sorted.bam \
									  FASTQ=/dev/stdout \
									  CLIPPING_ATTRIBUTE=XT \
									  CLIPPING_ACTION=N \
									  INTERLEAVE=true \
									  NON_PF=true \
									  TMP_DIR=$(TMPDIR) | \
									  bwa mem -M -t 12 -R \"@RG\tID:$$(*)\tLB:$$(*)\tPL:illumina\tSM:$$(*)\" \
									  -p $(REF_FASTA) /dev/stdin | \
									  $(JAVA) -Xmx8G -jar $(PICARD) SortSam \
									  I=/dev/stdin \
									  O=fgbio/$$(*).resorted.bam \
									  SORT_ORDER=coordinate \
									  TMP_DIR=$(TMPDIR)")
									  
fgbio/$1.fixed.bam : fgbio/%.resorted.bam
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
									  		   
fgbio/$1.intervals : fgbio/$1.fixed.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -w 1440,"set -o pipefail && \
									   							   		   /home/$(USER)/share/usr/jdk1.8.0_121/bin/java -Djava.io.tmpdir=$(TMPDIR) -Xms1G -Xmx12G -jar /home/$(USER)/share/usr/lib/java/GenomeAnalysisTK-3.7.jar \
									   							   		   -S LENIENT -T RealignerTargetCreator -I $$(^) -nt 8 -R $(REF_FASTA) -o $$(@) --known /home/$(USER)/share/reference/GATK_bundle/2.3/Mills_and_1000G_gold_standard.indels.b37.vcf.gz")

fgbio/$1.realn.bam : fgbio/$1.resorted.bam fgbio/$1.intervals
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -w 1440,"set -o pipefail && \
									   							   		   /home/$(USER)/share/usr/jdk1.8.0_121/bin/java -Djava.io.tmpdir=$(TMPDIR) -Xms1G -Xmx12G -jar /home/$(USER)/share/usr/lib/java/GenomeAnalysisTK-3.7.jar \
									   							   		   -S LENIENT -T IndelRealigner -I $$(<) -R $(REF_FASTA) -targetIntervals $$(<<) -o $$(@) --knownAlleles /home/brownd7/share/reference/GATK_bundle/2.3/Mills_and_1000G_gold_standard.indels.b37.vcf.gz")
									  		   
fgbio/$1.recal.grp : fgbio/$1.realn.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -w 1440,"set -o pipefail && \
																		   samtools index $$(<) && \
									   							   		   /home/$(USER)/share/usr/jdk1.8.0_121/bin/java -Djava.io.tmpdir=$(TMPDIR) -Xms1G -Xmx12G -jar /home/$(USER)/share/usr/lib/java/GenomeAnalysisTK-3.7.jar \
									   							   		   -S LENIENT -T BaseRecalibrator -R $(REF_FASTA) -knownSites /home/brownd7/share/reference/dbsnp_138.b37.gmaf.vcf.gz \
									   							   		   -I $$(<) -o $$(@)")

fgbio/$1.recal.bam : fgbio/$1.realn.bam fgbio/$1.recal.grp
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -w 1440,"set -o pipefail && \
									   							   		   /home/$(USER)/share/usr/jdk1.8.0_121/bin/java -Djava.io.tmpdir=$(TMPDIR) -Xms1G -Xmx12G -jar /home/$(USER)/share/usr/lib/java/GenomeAnalysisTK-3.7.jar \
									   							   		   -S LENIENT -T PrintReads -R $(REF_FASTA) -I $$(<) -BQSR $$(<<) -o $$(@)")

fgbio/$1.bam : fgbio/$1.recal.bam
	$$(call RUN, -c -n 1 -s 12G -m 18G -w 1440,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx16G -jar $$(PICARD_JAR) AddOrReplaceReadGroups \
												I=$$(<) \
												O=$$(@) \
												RGID=$1 \
												RGLB=$1 \
												RGPL=illumina \
												RGPU=NA \
												RGSM=$1 \
												TMP_DIR=$(TMPDIR) && \
												samtools index $$(@) && \
												cp fgbio/$1.bam.bai fgbio/$1.bai")

bam/%.bam : fgbio/%.bam
	$$(call RUN,-c -n 1 -s 4G -m 8G,"set -o pipefail && \
									 mkdir -p bam && \
									 cp fgbio/$$(*).bam bam/$$(*).bam && \
									 samtools index bam/$$(*).bam && \
									 cp bam/$$(*).bam.bai bam/$$(*).bai")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call bam-to-bam,$(sample))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
