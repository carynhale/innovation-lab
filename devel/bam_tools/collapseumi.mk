include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/collapse_umi.$(NOW)
PHONY += marianas

umi_collapse : $(foreach sample,$(SAMPLES),marianas/$(sample)/second-pass-alt-alleles.txt)

JAVA = $(HOME)/share/usr/jdk1.8.0_74/bin/java
MARIANAS_UMI_LENGTH ?= 3
MARIANAS_MIN_MAPQ ?= 1
MARIANAS_MIN_BAQ ?= 20
MARIANAS_MISMATCH ?= 0
MARIANAS_WOBBLE ?= 1
MARIANAS_MIN_CONSENSUS ?= 90
WALTZ_MIN_MAPQ ?= 20

define genotype-and-collapse
marianas/$1/$1.standard-pileup.txt : marianas/$1/$1.standard.bam
	$$(call RUN,-c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  cd marianas/$1 && \
									  $(JAVA) -server -Xms2G -Xmx8G -cp $(WALTZ) org.mskcc.juber.waltz.Waltz PileupMetrics $(WALTZ_MIN_MAPQ) $1.standard.bam $(REF_FASTA) $(WALTZ_BED_FILE) && \
									  cd ../..")
									  
marianas/$1/first-pass.mate-position-sorted.txt : marianas/$1/$1.standard-pileup.txt
	$$(call RUN,-c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  cd marianas/$1 && \
									  $(JAVA) -server -Xms2G -Xmx8G -cp $(MARIANAS) org.mskcc.marianas.umi.duplex.DuplexUMIBamToCollapsedFastqFirstPass \
									  $1.standard.bam \
									  $1.standard-pileup.txt \
									  $(MARIANAS_MIN_MAPQ) \
									  $(MARIANAS_MIN_BAQ) \
									  $(MARIANAS_MISMATCH) \
									  $(MARIANAS_WOBBLE) \
									  $(MARIANAS_MIN_CONSENSUS) \
									  $(REF_FASTA) && \
									  sort -n -s -S 6G -k 6 -k 8 first-pass.txt > first-pass.mate-position-sorted.txt && \
									  cd ../..")


marianas/$1/second-pass-alt-alleles.txt : marianas/$1/first-pass.mate-position-sorted.txt
	$$(call RUN,-c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  cd marianas/$1 && \
									  $(JAVA) -server -Xms2G -Xmx8G -cp $(MARIANAS) org.mskcc.marianas.umi.duplex.DuplexUMIBamToCollapsedFastqSecondPass \
									  $1.standard.bam \
									  $1.standard-pileup.txt \
									  $(MARIANAS_MIN_MAPQ) \
									  $(MARIANAS_MIN_BAQ) \
									  $(MARIANAS_MISMATCH) \
									  $(MARIANAS_WOBBLE) \
									  $(MARIANAS_MIN_CONSENSUS) \
									  $(REF_FASTA) \
									  first-pass.mate-position-sorted.txt && \
									  cd ../..")
									  									  
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call genotype-and-collapse,$(sample))))
