include innovation-lab/Makefile.inc
include innovation-lab/config/marianas.inc
include innovation-lab/config/waltz.inc
include innovation-lab/config/gatk.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/msk_access.$(NOW)

# MSK_ACCESS_WORKFLOW += copy_bam
# MSK_ACCESS_WORKFLOW += interval_metrics
# MSK_ACCESS_WORKFLOW += umi_qc
# MSK_ACCESS_WORKFLOW += plot_metrics
# MSK_ACCESS_WORKFLOW += cluster_samples

msk_access : $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)_R1.fastq.gz) \
		   	 $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)_R1_umi-clipped.fastq.gz) \
		   	 $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample).standard.bam) \
		   	 $(foreach sample,$(SAMPLES),marianas/$(sample)/second-pass-alt-alleles.txt) \
		   	 $(foreach sample,$(SAMPLES),marianas/$(sample)/timestamp)
		   	 
MARIANAS_UMI_LENGTH ?= 3
MARIANAS_MIN_MAPQ ?= 1
MARIANAS_MIN_BAQ ?= 20
MARIANAS_MISMATCH ?= 0
MARIANAS_WOBBLE ?= 1
MARIANAS_MIN_CONSENSUS ?= 90

WALTZ_MIN_MAPQ ?= 20
WALTZ_BED_FILE ?= $(HOME)/share/lib/bed_files/MSK-ACCESS-v1_0-probe-A.waltz.bed

BWA_ALN_OPTS ?= -M
BWAMEM_THREADS = 12
BWAMEM_MEM_PER_THREAD = 2G

SAMTOOLS_THREADS = 8
SAMTOOLS_MEM_THREAD = 2G

GATK_THREADS = 8
GATK_MEM_THREAD = 2G

define copy-fastq
marianas/$1/$1_R1.fastq.gz : $3
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
									 mkdir -p marianas/$1 && \
								     $(RSCRIPT) $(SCRIPTS_DIR)/fastq_tools/copy_fastq.R --sample_name $1 --fastq_files '$$^'")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call copy-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))
 
define clip-umi
marianas/$1/$1_R1_umi-clipped.fastq.gz : marianas/$1/$1_R1.fastq.gz
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  cd marianas/$1/ && \
									  $$(call MARIANAS_CMD,2G,8G) \
									  org.mskcc.marianas.umi.duplex.fastqprocessing.ProcessLoopUMIFastq \
									  $1_R1.fastq.gz $1_R2.fastq.gz \
									  $$(MARIANAS_UMI_LENGTH) && \
									  cd ../..")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call clip-umi,$(sample))))
	
define fastq-to-bam
marianas/$1/$1.bwamem.bam : marianas/$1/$1_R1_umi-clipped.fastq.gz
	$$(call RUN,-c -n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD),"set -o pipefail && \
																		   $$(BWA) mem -t $$(BWAMEM_THREADS) $$(BWA_ALN_OPTS) \
																		   -R \"@RG\tID:$1\tLB:$1\tPL:$$(SEQ_PLATFORM)\tSM:$1\" $$(REF_FASTA) marianas/$1/$1_R1_umi-clipped.fastq.gz marianas/$1/$1_R2_umi-clipped.fastq.gz | $$(SAMTOOLS) view -bhS - > $$(@)")
																		           
marianas/$1/$1.sorted.bam : marianas/$1/$1.bwamem.bam
	$$(call RUN,-c -n $(SAMTOOLS_THREADS) -s 1G -m $(SAMTOOLS_MEM_THREAD),"set -o pipefail && \
									  									   $$(SAMTOOLS) sort -@ $$(SAMTOOLS_THREADS) -m $$(SAMTOOLS_MEM_THREAD) $$(^) -o $$(@) -T $$(TMPDIR) && \
									  									   $$(SAMTOOLS) index $$(@) && \
									  									   cp marianas/$1/$1.sorted.bam.bai marianas/$1/$1.sorted.bai")
									  									   		   
marianas/$1/$1.fixed.bam : marianas/$1/$1.sorted.bam
	$$(call RUN,-c -n 1 -s 12G -m 16G,"set -o pipefail && \
									   $$(FIX_MATE) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   SORT_ORDER=coordinate \
									   COMPRESSION_LEVEL=0 \
									   CREATE_INDEX=true")

marianas/$1/$1.intervals : marianas/$1/$1.fixed.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD),"set -o pipefail && \
									   							   $$(call GATK_CMD,16G) \
									   							   -T RealignerTargetCreator \
									   							   -I $$(^) \
									   							   -nt $$(GATK_THREADS) \
									   							   -R $$(REF_FASTA) \
									   							   -o $$(@) \
									   							   -known $$(KNOWN_INDELS)")

marianas/$1/$1.realn.bam : marianas/$1/$1.sorted.bam marianas/$1/$1.intervals
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD),"set -o pipefail && \
									   							   $$(call GATK_CMD,16G) \
							   							   		   -T IndelRealigner \
							   							   		   -I $$(<) \
							   							   		   -R $$(REF_FASTA) \
							   							   		   -targetIntervals $$(<<) \
							   							   		   -o $$(@) \
									   							   -known $$(KNOWN_INDELS)")
									   							   
marianas/$1/$1.dedup.bam : marianas/$1/$1.realn.bam
	$$(call RUN, -c -n 1 -s 12G -m 18G,"set -o pipefail && \
										$$(MARK_DUP) \
										INPUT=$$(<) \
										OUTPUT=$$(@) \
										METRICS_FILE=marianas/$1/$1.txt \
										REMOVE_DUPLICATES=false \
										ASSUME_SORTED=true")
												
marianas/$1/$1.recal.grp : marianas/$1/$1.dedup.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD),"set -o pipefail && \
																   $$(SAMTOOLS) index $$(<) && \
									   							   $$(call GATK_CMD,16G) \
									   							   -T BaseRecalibrator \
									   							   -R $$(REF_FASTA) \
									   							   -knownSites $$(DBSNP) \
									   							   -I $$(<) \
									   							   -o $$(@)")

marianas/$1/$1.recal.bam : marianas/$1/$1.dedup.bam marianas/$1/$1.recal.grp
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD),"set -o pipefail && \
									   							   $$(call GATK_CMD,16G) \
									   							   -T PrintReads \
									   							   -R $$(REF_FASTA) \
									   							   -I $$(<) \
									   							   -BQSR $$(<<) \
									   							   -o $$(@)")

marianas/$1/$1.standard.bam : marianas/$1/$1.recal.bam
	$$(call RUN, -c -n 1 -s 12G -m 18G,"set -o pipefail && \
										$$(ADD_RG) \
										INPUT=$$(<) \
										OUTPUT=$$(@) \
										RGID=$1 \
										RGLB=$1 \
										RGPL=illumina \
										RGPU=NA \
										RGSM=$1 && \
										$$(SAMTOOLS) index $$(@) && \
										cp marianas/$1/$1.standard.bam.bai marianas/$1/$1.standard.bai")


endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fastq-to-bam,$(sample))))
		
define genotype-and-collapse
marianas/$1/$1.standard-pileup.txt : marianas/$1/$1.standard.bam
	$$(call RUN,-c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  cd marianas/$1 && \
									  $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1.standard.bam $$(REF_FASTA) $$(WALTZ_BED_FILE) && \
									  cd ../..")
									  
marianas/$1/first-pass.mate-position-sorted.txt : marianas/$1/$1.standard-pileup.txt
	$$(call RUN,-c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  cd marianas/$1 && \
									  $$(call MARIANAS_CMD,2G,8G) org.mskcc.marianas.umi.duplex.DuplexUMIBamToCollapsedFastqFirstPass \
									  $1.standard.bam \
									  $1.standard-pileup.txt \
									  $$(MARIANAS_MIN_MAPQ) \
									  $$(MARIANAS_MIN_BAQ) \
									  $$(MARIANAS_MISMATCH) \
									  $$(MARIANAS_WOBBLE) \
									  $$(MARIANAS_MIN_CONSENSUS) \
									  $$(REF_FASTA) && \
									  sort -n -s -S 6G -k 6 -k 8 first-pass.txt > first-pass.mate-position-sorted.txt && \
									  cd ../..")


marianas/$1/second-pass-alt-alleles.txt : marianas/$1/first-pass.mate-position-sorted.txt
	$$(call RUN,-c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  cd marianas/$1 && \
									  $$(call MARIANAS_CMD,2G,8G) org.mskcc.marianas.umi.duplex.DuplexUMIBamToCollapsedFastqSecondPass \
									  $1.standard.bam \
									  $1.standard-pileup.txt \
									  $$(MARIANAS_MIN_MAPQ) \
									  $$(MARIANAS_MIN_BAQ) \
									  $$(MARIANAS_MISMATCH) \
									  $$(MARIANAS_WOBBLE) \
									  $$(MARIANAS_MIN_CONSENSUS) \
									  $$(REF_FASTA) \
									  first-pass.mate-position-sorted.txt && \
									  cd ../..")
									  									  
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call genotype-and-collapse,$(sample))))
	
define fastq-to-collapsed-bam
marianas/$1/$1.collapsed.bwamem.bam : marianas/$1/second-pass-alt-alleles.txt
	$$(call RUN,-c -n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD),"set -o pipefail && \
																		   $$(BWA) mem -t $$(BWAMEM_THREADS) $$(BWA_ALN_OPTS) \
																		   -R \"@RG\tID:$1\tLB:$1\tPL:$$(SEQ_PLATFORM)\tSM:$1\" $$(REF_FASTA) marianas/$1/collapsed_R1_.fastq marianas/$1/collapsed_R2_.fastq | $$(SAMTOOLS) view -bhS - > $$(@)")
																		           
marianas/$1/$1.collapsed.sorted.bam : marianas/$1/$1.collapsed.bwamem.bam
	$$(call RUN,-c -n $(SAMTOOLS_THREADS) -s 1G -m $(SAMTOOLS_MEM_THREAD),"set -o pipefail && \
								  									   	   $$(SAMTOOLS) sort -@ $(SAMTOOLS_THREADS) -m $(SAMTOOLS_MEM_THREAD) $$(^) -o $$(@) -T $(TMPDIR) && \
								  									   	   $$(SAMTOOLS) index $$(@) && \
								  									   	   cp marianas/$1/$1.collapsed.sorted.bam.bai marianas/$1/$1.collapsed.sorted.bai")
								  									   		   
marianas/$1/$1.collapsed.fixed.bam : marianas/$1/$1.collapsed.sorted.bam
	$$(call RUN,-c -n 1 -s 12G -m 18G,"set -o pipefail && \
									   $$(FIX_MATE) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   SORT_ORDER=coordinate \
									   COMPRESSION_LEVEL=0 \
									   CREATE_INDEX=true")
									  		   
marianas/$1/$1.collapsed.intervals : marianas/$1/$1.collapsed.fixed.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD),"set -o pipefail && \
									   							   $$(call GATK_CMD,1G,12G)	\
									   							   -allowPotentiallyMisencodedQuals \
									   							   -T RealignerTargetCreator \
									   							   -I $$(^) \
									   							   -nt 8 \
									   							   -R $(REF_FASTA) \
									   							   -o $$(@) \
									   							    -known $$(KNOWN_INDELS)")

marianas/$1/$1.collapsed.realn.bam : marianas/$1/$1.collapsed.sorted.bam marianas/$1/$1.collapsed.intervals
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD),"set -o pipefail && \
									   							   $$(call GATK_CMD,1G,12G)	\
									   							   -allowPotentiallyMisencodedQuals \
									   							   -T IndelRealigner \
									   							   -I $$(<) \
									   							   -R $(REF_FASTA) \
									   							   -targetIntervals $$(<<) \
									   							   -o $$(@) \
									   							   -known $$(KNOWN_INDELS)")
									  		   
marianas/$1/$1.collapsed.bam : marianas/$1/$1.collapsed.realn.bam
	$$(call RUN, -c -n 1 -s 12G -m 18G,"set -o pipefail && \
										$$(ADD_RG) \
										INPUT=$$(<) \
										OUTPUT=$$(@) \
										RGID=$1 \
										RGLB=$1 \
										RGPL=illumina \
										RGPU=NA \
										RGSM=$1 && \
										$$(SAMTOOLS) index $$(@) && \
										cp marianas/$1/$1.collapsed.bam.bai marianas/$1/$1.collapsed.bai")
												
marianas/$1/timestamp : marianas/$1/$1.collapsed.bam
	$$(call RUN,-c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  cd marianas/$1 && \
									  $$(call MARIANAS_CMD,2G,8G) org.mskcc.marianas.umi.duplex.postprocessing.SeparateBams $1.collapsed.bam && \
									  echo 'Done!\n' > timestamp && \
									  cd ../..")
									  
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fastq-to-collapsed-bam,$(sample))))
	


# include modules/test/bam_tools/copybam.mk
# include modules/test/qc/intervalmetrics.mk
# include modules/test/qc/umiqc.mk
# include modules/test/qc/plotmetrics.mk
# include modules/test/qc/clustersamples.mk

..DUMMY := $(shell mkdir -p version; \
			 $(BWA) &> version/tmp.txt; \
			 head -3 version/tmp.txt | tail -2 > version/msk_access.txt; \
			 rm version/tmp.txt; \
			 $(SAMTOOLS) --version >> version/msk_access.txt; \
			 echo "gatk3" >> version/msk_access.txt; \
			 $(GATK) --version >> version/msk_access.txt; \
			 echo "picard" >> version/msk_access.txt; \
			 $(PICARD) MarkDuplicates --version &>> version/msk_access.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: msk_access
