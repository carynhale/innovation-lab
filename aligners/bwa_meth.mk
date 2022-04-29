include innovation-lab/Makefile.inc
include innovation-lab/config/gatk.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/bwa_meth.$(NOW)

bwa_meth : $(foreach sample,$(SAMPLES),bwameth/$(sample)/$(sample)_R1.fastq.gz) \
	   $(foreach sample,$(SAMPLES),bwameth/$(sample)/$(sample)_R2.fastq.gz) \
	   $(foreach sample,$(SAMPLES),bwameth/$(sample)/$(sample)_aln.bam) \
	   $(foreach sample,$(SAMPLES),bwameth/$(sample)/$(sample)_aln_srt.bam) \
	   $(foreach sample,$(SAMPLES),bwameth/$(sample)/$(sample)_aln_srt_MD.bam)
#	   $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD.intervals) \
#	   $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR.bam) \
#	   $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX.bam) \
#	   $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX_BR.bam) \
#	   $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp.bam) \
#	   $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC.bam) \
#	   $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC.duplex_umi_counts.txt) \
#	   $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA.bam) \
#	   $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG.bam) \
#	   $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG.intervals) \
#	   $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR.bam) \
#	   $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam) \
#	   $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bam) \
#	   $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bam) \
#	   $(foreach sample,$(SAMPLES),bam/$(sample)_cl_aln_srt_MD_IR_FX_BR.bam)

BWAMETH_GENOME = $(REF_DIR)/IDT_oligo/bwa_meth/idt_oligo.fasta

BWAMETH_THREADS = 12
BWAMETH_MEM_PER_THREAD = 2G

SAMTOOLS_THREADS = 8
SAMTOOLS_MEM_THREAD = 2G

GATK_THREADS = 8
GATK_MEM_THREAD = 2G

define merge-fastq
bwameth/$1/$1_R1.fastq.gz : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G,"zcat $$(^) | gzip -c > $$(@)")
	
bwameth/$1/$1_R2.fastq.gz : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G,"zcat $$(^) | gzip -c > $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merge-fastq,$(sample),$(split.$(sample)))))

define fastq-2-bam
bwameth/$1/$1_aln.bam : bwameth/$1/$1_R1.fastq.gz bwameth/$1/$1_R2.fastq.gz
	$$(call RUN,-c -n $(BWAMETH_THREADS) -s 1G -m $(BWAMETH_MEM_PER_THREAD) -v $(BWAMETH_ENV),"set -o pipefail && \
												   $$(BWAMETH) \
												   --threads $$(BWAMETH_THREADS) \
												   --reference $$(BWAMETH_GENOME) \
												   $$(<) $$(<<) | \
												   $$(SAMTOOLS) view -bhS - > $$@")

bwameth/$1/$1_aln_srt.bam : bwameth/$1/$1_aln.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
					  $$(SORT_SAM) \
					  INPUT=$$(<) \
					  OUTPUT=$$(@) \
					  SORT_ORDER=coordinate \
					  $$(SAMTOOLS) \
					  index \
					  $$(@) && \
					  cp bwameth/$1/$1_aln_srt.bam.bai bwameth/$1/$1_aln_srt.bai")
					 
bwameth/$1/$1_aln_srt_MD.bam : bwameth/$1/$1_aln_srt.bam
	$$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
					    $$(MARK_DUP) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@) \
					    METRICS_FILE=bwameth/$1/$1_aln_srt.txt \
					    REMOVE_DUPLICATES=false \
					    ASSUME_SORTED=true && \
					    $$(SAMTOOLS) index $$(@) && \
					    cp bwameth/$1/$1_aln_srt_MD.bam.bai bwameth/$1/$1_aln_srt_MD.bai")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call fastq-2-bam,$(sample))))
	


define merge-bams


fgbio/$1/$1_cl.fastq.gz : fgbio/$1/$1_fq_srt.bam
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



fgbio/$1/$1_cl_aln_srt_MD.intervals : fgbio/$1/$1_cl_aln_srt_MD.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -v $(GATK_ENV),"set -o pipefail && \
										      $$(call GATK_CMD,16G) \
										      -T RealignerTargetCreator \
										      -I $$(^) \
										      -nt $$(GATK_THREADS) \
										      -R $$(REF_FASTA) \
										      -o $$(@) \
										      -known $$(KNOWN_INDELS)")

fgbio/$1/$1_cl_aln_srt_MD_IR.bam : fgbio/$1/$1_cl_aln_srt_MD.bam fgbio/$1/$1_cl_aln_srt_MD.intervals
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -v $(GATK_ENV),"set -o pipefail && \
										      $$(call GATK_CMD,16G) \
										      -T IndelRealigner \
										      -I $$(<) \
										      -R $$(REF_FASTA) \
										      -targetIntervals $$(<<) \
										      -o $$(@) \
										      -known $$(KNOWN_INDELS)")
									   							   
fgbio/$1/$1_cl_aln_srt_MD_IR_FX.bam : fgbio/$1/$1_cl_aln_srt_MD_IR.bam
	$$(call RUN,-c -n 1 -s 24G -m 36G,"set -o pipefail && \
					   $$(FIX_MATE) \
					   INPUT=$$(<) \
					   OUTPUT=$$(@) \
					   SORT_ORDER=coordinate \
					   COMPRESSION_LEVEL=0 \
					   CREATE_INDEX=true")

fgbio/$1/$1_cl_aln_srt_MD_IR_FX.grp : fgbio/$1/$1_cl_aln_srt_MD_IR_FX.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -v $(GATK_ENV),"set -o pipefail && \
										      $$(SAMTOOLS) index $$(<) && \
										      $$(call GATK_CMD,16G) \
										      -T BaseRecalibrator \
										      -R $$(REF_FASTA) \
										      -knownSites $$(DBSNP) \
										      -I $$(<) \
										      -o $$(@)")

fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX.bam fgbio/$1/$1_cl_aln_srt_MD_IR_FX.grp
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -v $(GATK_ENV),"set -o pipefail && \
										      $$(call GATK_CMD,16G) \
										      -T PrintReads \
										      -R $$(REF_FASTA) \
										      -I $$(<) \
										      -BQSR $$(<<) \
										      -o $$(@)")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call merge-bams,$(sample))))
	
define create-consensus
fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
					  $$(call FGBIO_CMD,2G,8G) \
					  GroupReadsByUmi \
					  --strategy=paired \
					  --input=$$(<) \
					  --output=$$(@) \
					  --family-size-histogram=fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp.txt")
									  
fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD),"set -o pipefail && \
								       $$(call FGBIO_CMD,2G,16G) \
								       CallDuplexConsensusReads \
								       --min-reads 1 1 0 \
								       --input $$(<) \
								       --output $$(@) \
								       --threads $$(GATK_THREADS)")

fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC.duplex_umi_counts.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
					  $$(call FGBIO_CMD,2G,16G) \
					  CollectDuplexSeqMetrics \
					  --input $$(<) \
					  --output fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC \
					  --duplex-umi-counts true \
					  --intervals $$(TARGETS_LIST)")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call create-consensus,$(sample))))
	
define align-consensus
fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC.bam
	$$(call RUN,-c -n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD),"set -o pipefail && \
									       $$(SAM_TO_FASTQ) \
									       INPUT=$$(<) \
									       FASTQ=/dev/stdout \
									       INTERLEAVE=true | \
									       $$(BWA) mem -p -t $$(BWAMEM_THREADS) $$(REF_FASTA) /dev/stdin | \
									       $$(MERGE_ALIGNMENTS) \
									       UNMAPPED=$$(<) \
									       ALIGNED=/dev/stdin \
									       OUTPUT=$$(@) \
									       REFERENCE_SEQUENCE=$$(REF_FASTA) \
									       MAX_GAPS=-1 \
									       ORIENTATIONS=FR")

fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(ADD_RG) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@) \
					    RGID=$1 \
					    RGLB=$1 \
					    RGPL=illumina \
					    RGPU=NA \
					    RGSM=$1 \
					    SORT_ORDER=coordinate \
					    COMPRESSION_LEVEL=0 && \
					    $$(SAMTOOLS) index $$(@) && \
					    cp fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG.bam.bai fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG.bai")

fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG.intervals : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -v $(GATK_ENV),"set -o pipefail && \
										      $$(call GATK_CMD,16G) \
										      -T RealignerTargetCreator \
										      -I $$(^) \
										      -nt $$(GATK_THREADS) \
										      -R $$(REF_FASTA) \
										      -o $$(@) \
										      -known $$(KNOWN_INDELS) \
										      --allow_potentially_misencoded_quality_scores")

fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG.bam fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG.intervals
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -v $(GATK_ENV),"set -o pipefail && \
										      $$(call GATK_CMD,16G) \
										      -T IndelRealigner \
										      -I $$(<) \
										      -R $$(REF_FASTA) \
										      -targetIntervals $$(<<) \
										      -o $$(@) \
										      -known $$(KNOWN_INDELS) \
										      --allow_potentially_misencoded_quality_scores")
									   							   
fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR.bam
	$$(call RUN,-c -n 1 -s 24G -m 36G,"set -o pipefail && \
					   $$(FIX_MATE) \
					   INPUT=$$(<) \
					   OUTPUT=$$(@) \
					   SORT_ORDER=coordinate \
					   COMPRESSION_LEVEL=0 \
					   CREATE_INDEX=true")
									   
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call align-consensus,$(sample))))
	
define filter-consensus
fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
					  $$(PYTHON) $$(SCRIPTS_DIR)/bam_tools/create_simplex_bam_from_consensus.py \
					  $$(<) \
					  $$(@)")

fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
					  $$(call FGBIO_CMD,2G,16G) \
					  FilterConsensusReads \
					  --input $$(<) \
					  --output $$(@) \
					  --ref $$(REF_FASTA) \
					  --min-reads=2 1 1 \
					  --min-base-quality=30 \
					  --reverse-per-base-tags=true")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call filter-consensus,$(sample))))
	
define copy-bam
bam/$1_cl_aln_srt_MD_IR_FX_BR.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR.bam
	$$(call RUN, -c -s 2G -m 4G,"set -o pipefail && \
				     cp $$(<) $$(@) && \
				     $$(SAMTOOLS) index $$(@) && \
				     cp bam/$1_cl_aln_srt_MD_IR_FX_BR.bam.bai bam/$1_cl_aln_srt_MD_IR_FX_BR.bai")

bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam
	$$(call RUN, -c -s 2G -m 4G,"set -o pipefail && \
				     cp $$(<) $$(@) && \
				     $$(SAMTOOLS) index $$(@) && \
				     cp bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam.bai bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bai")
												
bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bam
	$$(call RUN, -c -s 2G -m 4G,"set -o pipefail && \
				     cp $$(<) $$(@) && \
				     $$(SAMTOOLS) index $$(@) && \
				     cp bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bam.bai bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bai")
												
bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bam
	$$(call RUN, -c -s 2G -m 4G,"set -o pipefail && \
				     cp $$(<) $$(@) && \
				     $$(SAMTOOLS) index $$(@) && \
				     cp bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bam.bai bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bai")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call copy-bam,$(sample))))
	
define picard-metrics-standard
metrics/$1_cl_aln_srt_MD_IR_FX_BR.idx_stats.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(BAM_INDEX) \
					    INPUT=$$(<) \
					    > $$(@)")
									   
metrics/$1_cl_aln_srt_MD_IR_FX_BR.aln_metrics.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(COLLECT_ALIGNMENT_METRICS) \
					    REFERENCE_SEQUENCE=$$(REF_FASTA) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@)")
									   
metrics/$1_cl_aln_srt_MD_IR_FX_BR.insert_metrics.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(COLLECT_INSERT_METRICS) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@) \
					    HISTOGRAM_FILE=metrics/$1_cl_aln_srt_MD_IR_FX_BR.insert_metrics.pdf \
					    MINIMUM_PCT=0.5")
									   
metrics/$1_cl_aln_srt_MD_IR_FX_BR.oxog_metrics.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(COLLECT_OXOG_METRICS) \
					    REFERENCE_SEQUENCE=$$(REF_FASTA) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@)")
					    
metrics/$1_cl_aln_srt_MD_IR_FX_BR.hs_metrics.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					   $$(CALC_HS_METRICS) \
					   REFERENCE_SEQUENCE=$$(REF_FASTA) \
					   INPUT=$$(<) \
					   OUTPUT=$$(@) \
					   BAIT_INTERVALS=$$(TARGETS_LIST) \
					   TARGET_INTERVALS=$$(TARGETS_LIST) \
					   MINIMUM_MAPPING_QUALITY=0 \
					   PER_TARGET_COVERAGE=metrics/$1_cl_aln_srt_MD_IR_FX_BR.hs_metrics_target.txt \
					   PER_BASE_COVERAGE=metrics/$1_cl_aln_srt_MD_IR_FX_BR.hs_metrics_base.txt")
					   
metrics/$1_cl_aln_srt_MD_IR_FX_BR.gc_metrics_summary.txt : bam/$1_cl_aln_srt_MD_IR_FX_BR.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
					   $$(COLLECT_GC_BIAS) \
					   INPUT=$$(<) \
					   OUTPUT=metrics/$1_cl_aln_srt_MD_IR_FX_BR.gc_metrics.txt \
					   CHART_OUTPUT=metrics/$1_cl_aln_srt_MD_IR_FX_BR.gc_metrics.pdf \
					   REFERENCE_SEQUENCE=$$(REF_FASTA) \
					   SUMMARY_OUTPUT=$$(@)")

endef
$(foreach sample,$(SAMPLES),\
 		$(eval $(call picard-metrics-standard,$(sample))))
		
define picard-metrics-collapsed
metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.idx_stats.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(BAM_INDEX) \
					    INPUT=$$(<) \
					    > $$(@)")
									   
metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.aln_metrics.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(COLLECT_ALIGNMENT_METRICS) \
					    REFERENCE_SEQUENCE=$$(REF_FASTA) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@)")
									   
metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.insert_metrics.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(COLLECT_INSERT_METRICS) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@) \
					    HISTOGRAM_FILE=metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.insert_metrics.pdf \
					    MINIMUM_PCT=0.5")
									   
metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.oxog_metrics.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(COLLECT_OXOG_METRICS) \
					    REFERENCE_SEQUENCE=$$(REF_FASTA) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@)")
					    
metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.hs_metrics.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					   $$(CALC_HS_METRICS) \
					   REFERENCE_SEQUENCE=$$(REF_FASTA) \
					   INPUT=$$(<) \
					   OUTPUT=$$(@) \
					   BAIT_INTERVALS=$$(TARGETS_LIST) \
					   TARGET_INTERVALS=$$(TARGETS_LIST) \
					   MINIMUM_MAPPING_QUALITY=0 \
					   PER_TARGET_COVERAGE=metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.hs_metrics_target.txt \
					   PER_BASE_COVERAGE=metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.hs_metrics_base.txt")

metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.gc_metrics_summary.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
					   $$(COLLECT_GC_BIAS) \
					   INPUT=$$(<) \
					   OUTPUT=metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.gc_metrics.txt \
					   CHART_OUTPUT=metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.gc_metrics.pdf \
					   REFERENCE_SEQUENCE=$$(REF_FASTA) \
					   SUMMARY_OUTPUT=$$(@)")

endef
$(foreach sample,$(SAMPLES),\
 		$(eval $(call picard-metrics-collapsed,$(sample))))
					    
define picard-metrics-simplex
metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.idx_stats.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(BAM_INDEX) \
					    INPUT=$$(<) \
					    > $$(@)")
									   
metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.aln_metrics.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(COLLECT_ALIGNMENT_METRICS) \
					    REFERENCE_SEQUENCE=$$(REF_FASTA) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@)")
									   
metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.insert_metrics.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(COLLECT_INSERT_METRICS) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@) \
					    HISTOGRAM_FILE=metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.insert_metrics.pdf \
					    MINIMUM_PCT=0.5")
									   
metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.oxog_metrics.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(COLLECT_OXOG_METRICS) \
					    REFERENCE_SEQUENCE=$$(REF_FASTA) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@)")
					    
metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.hs_metrics.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					   $$(CALC_HS_METRICS) \
					   REFERENCE_SEQUENCE=$$(REF_FASTA) \
					   INPUT=$$(<) \
					   OUTPUT=$$(@) \
					   BAIT_INTERVALS=$$(TARGETS_LIST) \
					   TARGET_INTERVALS=$$(TARGETS_LIST) \
					   MINIMUM_MAPPING_QUALITY=0 \
					   PER_TARGET_COVERAGE=metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.hs_metrics_target.txt \
					   PER_BASE_COVERAGE=metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.hs_metrics_base.txt")

metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.gc_metrics_summary.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
					   $$(COLLECT_GC_BIAS) \
					   INPUT=$$(<) \
					   OUTPUT=metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.gc_metrics.txt \
					   CHART_OUTPUT=metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.gc_metrics.pdf \
					   REFERENCE_SEQUENCE=$$(REF_FASTA) \
					   SUMMARY_OUTPUT=$$(@)")

endef
$(foreach sample,$(SAMPLES),\
 		$(eval $(call picard-metrics-simplex,$(sample))))
		
define picard-metrics-duplex
metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.idx_stats.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(BAM_INDEX) \
					    INPUT=$$(<) \
					    > $$(@)")
									   
metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.aln_metrics.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(COLLECT_ALIGNMENT_METRICS) \
					    REFERENCE_SEQUENCE=$$(REF_FASTA) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@)")
									   
metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.insert_metrics.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(COLLECT_INSERT_METRICS) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@) \
					    HISTOGRAM_FILE=metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.insert_metrics.pdf \
					    MINIMUM_PCT=0.5")
									   
metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.oxog_metrics.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					    $$(COLLECT_OXOG_METRICS) \
					    REFERENCE_SEQUENCE=$$(REF_FASTA) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@)")
					    
metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.hs_metrics.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G,"set -o pipefail && \
					   $$(CALC_HS_METRICS) \
					   REFERENCE_SEQUENCE=$$(REF_FASTA) \
					   INPUT=$$(<) \
					   OUTPUT=$$(@) \
					   BAIT_INTERVALS=$$(TARGETS_LIST) \
					   TARGET_INTERVALS=$$(TARGETS_LIST) \
					   MINIMUM_MAPPING_QUALITY=0 \
					   PER_TARGET_COVERAGE=metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.hs_metrics_target.txt \
					   PER_BASE_COVERAGE=metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.hs_metrics_base.txt")

metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.gc_metrics_summary.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
					   $$(COLLECT_GC_BIAS) \
					   INPUT=$$(<) \
					   OUTPUT=metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.gc_metrics.txt \
					   CHART_OUTPUT=metrics/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.gc_metrics.pdf \
					   REFERENCE_SEQUENCE=$$(REF_FASTA) \
					   SUMMARY_OUTPUT=$$(@)")

endef
$(foreach sample,$(SAMPLES),\
 		$(eval $(call picard-metrics-duplex,$(sample))))
		
summary/idx_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR.idx_stats.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.idx_stats.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.idx_stats.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/fgbio_access.R --option 1 --sample_names '$(SAMPLES)'")

summary/aln_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR.aln_metrics.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.aln_metrics.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.aln_metrics.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.aln_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/fgbio_access.R --option 2 --sample_names '$(SAMPLES)'")

summary/insert_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR.insert_metrics.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.insert_metrics.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.insert_metrics.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/fgbio_access.R --option 3 --sample_names '$(SAMPLES)'")

summary/oxog_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR.oxog_metrics.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.oxog_metrics.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.oxog_metrics.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.oxog_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/fgbio_access.R --option 4 --sample_names '$(SAMPLES)'")

summary/hs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR.hs_metrics.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.hs_metrics.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.hs_metrics.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.hs_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/fgbio_access.R --option 5 --sample_names '$(SAMPLES)'")
					  
summary/umi_counts.txt : $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC.duplex_umi_counts.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/fgbio_access.R --option 6 --sample_names '$(SAMPLES)'")
									  
summary/umi_duplex_counts.txt : $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC.duplex_umi_counts.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/fgbio_access.R --option 7 --sample_names '$(SAMPLES)'")
									  
summary/all_family_sizes.txt : $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC.duplex_umi_counts.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/fgbio_access.R --option 8 --sample_names '$(SAMPLES)'")
									  
summary/duplex_family_sizes.txt : $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC.duplex_umi_counts.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/fgbio_access.R --option 9 --sample_names '$(SAMPLES)'")

summary/duplex_yield_metrics.txt : $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC.duplex_umi_counts.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/fgbio_access.R --option 10 --sample_names '$(SAMPLES)'")

summary/gc_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR.gc_metrics_summary.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.gc_metrics_summary.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.gc_metrics_summary.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.gc_metrics_summary.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/fgbio_access.R --option 11 --sample_names '$(SAMPLES)'")


..DUMMY := $(shell mkdir -p version; \
	     $(JAVA8) -jar $(FGBIO) --help &> version/bwa_meth.txt; \
	     echo "picard" >> version/bwa_meth.txt; \
	     $(PICARD) SortSam --version &>> version/bwa_meth.txt; \
	     $(PICARD) MarkIlluminaAdapters --version &>> version/bwa_meth.txt; \
	     $(PICARD) SamToFastq --version &>> version/bwa_meth.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: bwa_meth
