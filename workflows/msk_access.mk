include innovation-lab/Makefile.inc
include innovation-lab/config/marianas.inc
include innovation-lab/config/waltz.inc
include innovation-lab/config/gatk.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/msk_access.$(NOW)

# MSK_ACCESS_WORKFLOW += plot_metrics
# MSK_ACCESS_WORKFLOW += cluster_samples
# include modules/test/qc/plotmetrics.mk
# include modules/test/qc/clustersamples.mk

msk_access : $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)_R1.fastq.gz) \
		   	 $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)_R1_umi-clipped.fastq.gz) \
		   	 $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample).standard.bam) \
		   	 $(foreach sample,$(SAMPLES),marianas/$(sample)/second-pass-alt-alleles.txt) \
		   	 $(foreach sample,$(SAMPLES),marianas/$(sample)/timestamp) \
			 $(foreach sample,$(SAMPLES),bam/$(sample)-standard.bam) \
			 $(foreach sample,$(SAMPLES),bam/$(sample)-unfiltered.bam) \
			 $(foreach sample,$(SAMPLES),bam/$(sample)-simplex.bam) \
			 $(foreach sample,$(SAMPLES),bam/$(sample)-duplex.bam) \
			 $(foreach sample,$(SAMPLES),marianas/$(sample)/family-sizes.txt) \
			 metrics/summary/umi_frequencies.tsv \
			 metrics/summary/umi_composite.tsv \
			 metrics/summary/umi_families.tsv \
 			 $(foreach sample,$(SAMPLES),metrics/standard/$(sample)-pileup.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/simplex/$(sample)-pileup.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/duplex/$(sample)-pileup.txt) \
			 $(foreach sample,$(SAMPLES),metrics/standard/$(sample).A.ontarget.txt) \
			 $(foreach sample,$(SAMPLES),metrics/standard/$(sample).B.ontarget.txt) \
			 $(foreach sample,$(SAMPLES),metrics/standard/$(sample).AB.offtarget.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/standard/$(sample).idx_stats.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/standard/$(sample).aln_metrics.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/standard/$(sample).insert_metrics.txt) \
  			 $(foreach sample,$(SAMPLES),metrics/standard/$(sample).oxog_metrics.txt) \
  			 $(foreach sample,$(SAMPLES),metrics/standard/$(sample).probe-A.hs_metrics.txt) \
  			 $(foreach sample,$(SAMPLES),metrics/standard/$(sample).probe-B.hs_metrics.txt) \
  			 $(foreach sample,$(SAMPLES),metrics/standard/$(sample).probe-A.hs_metrics-nodedup.txt) \
  			 $(foreach sample,$(SAMPLES),metrics/standard/$(sample).probe-B.hs_metrics-nodedup.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample).idx_stats.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample).aln_metrics.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample).insert_metrics.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample).probe-A.hs_metrics.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample).probe-B.hs_metrics.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/simplex/$(sample).idx_stats.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/simplex/$(sample).aln_metrics.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/simplex/$(sample).insert_metrics.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/simplex/$(sample).probe-A.hs_metrics.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/simplex/$(sample).probe-B.hs_metrics.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/duplex/$(sample).idx_stats.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/duplex/$(sample).aln_metrics.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/duplex/$(sample).insert_metrics.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/duplex/$(sample).probe-A.hs_metrics.txt) \
 			 $(foreach sample,$(SAMPLES),metrics/duplex/$(sample).probe-B.hs_metrics.txt) \
 			 metrics/standard/metrics_idx.tsv \
 			 metrics/unfiltered/metrics_idx.tsv \
 			 metrics/simplex/metrics_idx.tsv \
 			 metrics/duplex/metrics_idx.tsv \
 			 metrics/standard/metrics_aln.tsv \
 			 metrics/unfiltered/metrics_aln.tsv \
 			 metrics/simplex/metrics_aln.tsv \
 			 metrics/duplex/metrics_aln.tsv \
 			 metrics/standard/metrics_insert.tsv \
 			 metrics/unfiltered/metrics_insert.tsv \
 			 metrics/simplex/metrics_insert.tsv \
 			 metrics/duplex/metrics_insert.tsv \
 			 metrics/standard/metrics_insert_distribution.tsv \
 			 metrics/unfiltered/metrics_insert_distribution.tsv \
 			 metrics/simplex/metrics_insert_distribution.tsv \
 			 metrics/duplex/metrics_insert_distribution.tsv \
 			 metrics/standard/metrics_hs.tsv \
 			 metrics/unfiltered/metrics_hs.tsv \
 			 metrics/simplex/metrics_hs.tsv \
 			 metrics/duplex/metrics_hs.tsv \
 			 metrics/standard/metrics_oxog.tsv \
 			 metrics/summary/metrics_idx.tsv \
 			 metrics/summary/metrics_aln.tsv \
 			 metrics/summary/metrics_insert.tsv \
 			 metrics/summary/metrics_insert_distribution.tsv \
 			 metrics/summary/metrics_hs.tsv \
 			 metrics/summary/metrics_ts.tsv

WALTZ_BED_FILE ?= $(HOME)/share/lib/bed_files/MSK-ACCESS-v1_0-probe-A.sorted.bed
UMI_QC_BED_FILE_A ?= $(HOME)/share/lib/bed_files/MSK-ACCESS-v1_0-probe-A.sorted.bed
UMI_QC_BED_FILE_B ?= $(HOME)/share/lib/bed_files/MSK-ACCESS-v1_0-probe-B.sorted.bed
OFF_TARGET_FILE_AB ?= $(HOME)/share/lib/bed_files/MSK-ACCESS-v1_0-probe-AB.offtarget.bed
POOL_A_TARGET_FILE ?= $(HOME)/share/lib/resource_files/MSK-ACCESS-v1_0-probe-A.sorted.list
POOL_B_TARGET_FILE ?= $(HOME)/share/lib/resource_files/MSK-ACCESS-v1_0-probe-B.sorted.list

MARIANAS_UMI_LENGTH ?= 3
MARIANAS_MIN_MAPQ ?= 1
MARIANAS_MIN_BAQ ?= 20
MARIANAS_MISMATCH ?= 0
MARIANAS_WOBBLE ?= 1
MARIANAS_MIN_CONSENSUS ?= 90

WALTZ_MIN_MAPQ ?= 20

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
									  cut -f 1,2,3 $$(WALTZ_BED_FILE) > .bed && \
									  $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1.standard.bam $$(REF_FASTA) .bed && \
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
								 cp marianas/$1/$1.collapsed-simplex.bam $$(@) && \
								 cp marianas/$1/$1.collapsed-simplex.bai bam/$1-simplex.bam.bai && \
								 cp marianas/$1/$1.collapsed-simplex.bai bam/$1-simplex.bai")
												
bam/$1-duplex.bam : marianas/$1/timestamp
	$$(call RUN, -c -s 2G -m 4G,"set -o pipefail && \
								 cp marianas/$1/$1.collapsed-duplex.bam $$(@) && \
								 cp marianas/$1/$1.collapsed-duplex.bai bam/$1-duplex.bam.bai && \
								 cp marianas/$1/$1.collapsed-duplex.bai bam/$1-duplex.bai")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call copy-to-bam,$(sample))))

define family-size-metric
marianas/$1/family-sizes.txt : marianas/$1/second-pass-alt-alleles.txt
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   cd marianas/$1 && \
									   source ../../$$(SCRIPTS_DIR)/qc/umi_metrics.sh $$(UMI_QC_BED_FILE_A) $$(UMI_QC_BED_FILE_B) $1")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call family-size-metric,$(sample))))

define pileup-metric
metrics/standard/$1-pileup.txt : bam/$1-standard.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
								 cd metrics/standard && \
								 ln -sf ../../bam/$1-standard.bam $1.bam && \
								 ln -sf ../../bam/$1-standard.bam.bai $1.bam.bai && \
								 ln -sf ../../bam/$1-standard.bai $1.bai && \
								 if [[ ! -f '.bed' ]]; then cut -f 4 $$(WALTZ_BED_FILE) | paste -d '\t' $$(WALTZ_BED_FILE) - > .bed; fi && \
								 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1.bam $$(REF_FASTA) .bed && \
								 unlink $1.bam && \
								 unlink $1.bam.bai && \
								 unlink $1.bai && \
								 rm -rf .bed && \
								 cd ../..")
									 
metrics/simplex/$1-pileup.txt : bam/$1-simplex.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
								 cd metrics/simplex && \
								 ln -sf ../../bam/$1-simplex.bam $1.bam && \
								 ln -sf ../../bam/$1-simplex.bam.bai $1.bam.bai && \
								 ln -sf ../../bam/$1-simplex.bai $1.bai && \
								 if [[ ! -f '.bed' ]]; then cut -f 4 $$(WALTZ_BED_FILE) | paste -d '\t' $$(WALTZ_BED_FILE) - > .bed; fi && \
								 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1.bam $$(REF_FASTA) .bed && \
								 unlink $1.bam && \
								 unlink $1.bam.bai && \
								 unlink $1.bai && \
								 rm -rf .bed && \
								 cd ../..")
								 
metrics/duplex/$1-pileup.txt : bam/$1-duplex.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
								 cd metrics/duplex && \
								 ln -sf ../../bam/$1-duplex.bam $1.bam && \
								 ln -sf ../../bam/$1-duplex.bam.bai $1.bam.bai && \
								 ln -sf ../../bam/$1-duplex.bai $1.bai && \
								 if [[ ! -f '.bed' ]]; then cut -f 4 $$(WALTZ_BED_FILE) | paste -d '\t' $$(WALTZ_BED_FILE) - > .bed; fi && \
								 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1.bam $$(REF_FASTA) .bed && \
								 unlink $1.bam && \
								 unlink $1.bam.bai && \
								 unlink $1.bai && \
								 rm -rf .bed && \
								 cd ../..")
								 
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call pileup-metric,$(sample))))
				   
define coverage-metric
metrics/standard/$1.A.ontarget.txt : bam/$1-standard.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
								 $$(SAMTOOLS) view -L $$(UMI_QC_BED_FILE_A) $$(<) -b > metrics/standard/$1-ontarget-A.bam && \
								 $$(SAMTOOLS) index metrics/standard/$1-ontarget-A.bam && \
								 $$(BAM_INDEX) \
								 INPUT=metrics/standard/$1-ontarget-A.bam \
								 > $$(@) && \
								 rm -rf metrics/standard/$1-ontarget-A.bam && \
								 rm -rf metrics/standard/$1-ontarget-A.bam.bai")
									 
metrics/standard/$1.B.ontarget.txt : bam/$1-standard.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
								 $$(SAMTOOLS) view -L $$(UMI_QC_BED_FILE_B) $$(<) -b > metrics/standard/$1-ontarget-B.bam && \
								 $$(SAMTOOLS) index metrics/standard/$1-ontarget-B.bam && \
								 $$(BAM_INDEX) \
								 INPUT=metrics/standard/$1-ontarget-B.bam \
								 > $$(@) && \
								 rm -rf metrics/standard/$1-ontarget-B.bam && \
								 rm -rf metrics/standard/$1-ontarget-B.bam.bai")
	
metrics/standard/$1.AB.offtarget.txt : bam/$1-standard.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
								 $$(SAMTOOLS) view -L $$(OFF_TARGET_FILE_AB) $$(<) -b > metrics/standard/$1-offtarget-AB.bam && \
								 $$(SAMTOOLS) index metrics/standard/$1-offtarget-AB.bam && \
								 $$(BAM_INDEX) \
								 INPUT=metrics/standard/$1-offtarget-AB.bam \
								 > $$(@) && \
								 rm -rf metrics/standard/$1-offtarget-AB.bam && \
								 rm -rf metrics/standard/$1-offtarget-AB.bam.bai")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call coverage-metric,$(sample))))
 
define picard-metrics-standard
metrics/standard/$1.idx_stats.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(BAM_INDEX) \
									   INPUT=$$(<) \
									   > $$(@)")
									   
metrics/standard/$1.aln_metrics.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(COLLECT_ALIGNMENT_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@)")
									   
metrics/standard/$1.insert_metrics.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(COLLECT_INSERT_METRICS) \
 									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   HISTOGRAM_FILE=metrics/standard/$1.insert_metrics.pdf \
									   MINIMUM_PCT=0.5")
									   
metrics/standard/$1.oxog_metrics.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(COLLECT_OXOG_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@)")

metrics/standard/$1.probe-A.hs_metrics.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(CALC_HS_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   BAIT_INTERVALS=$$(POOL_A_TARGET_FILE) \
									   TARGET_INTERVALS=$$(POOL_A_TARGET_FILE)")
 												
metrics/standard/$1.probe-B.hs_metrics.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(CALC_HS_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   BAIT_INTERVALS=$$(POOL_B_TARGET_FILE) \
									   TARGET_INTERVALS=$$(POOL_B_TARGET_FILE)")
									   
metrics/standard/$1.probe-A.hs_metrics-nodedup.txt : marianas/$1/$1.realn.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(CALC_HS_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   BAIT_INTERVALS=$$(POOL_A_TARGET_FILE) \
									   TARGET_INTERVALS=$$(POOL_A_TARGET_FILE)")
												
metrics/standard/$1.probe-B.hs_metrics-nodedup.txt : marianas/$1/$1.realn.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(CALC_HS_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   BAIT_INTERVALS=$$(POOL_B_TARGET_FILE) \
									   TARGET_INTERVALS=$$(POOL_B_TARGET_FILE)")

endef
$(foreach sample,$(SAMPLES),\
 		$(eval $(call picard-metrics-standard,$(sample))))
 		
define picard-metrics-unfiltered
metrics/unfiltered/$1.idx_stats.txt : bam/$1-unfiltered.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(BAM_INDEX) \
									   INPUT=$$(<) \
									   > $$(@)")
									   
metrics/unfiltered/$1.aln_metrics.txt : bam/$1-unfiltered.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(COLLECT_ALIGNMENT_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@)")
									   
metrics/unfiltered/$1.insert_metrics.txt : bam/$1-unfiltered.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(COLLECT_INSERT_METRICS) \
 									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   HISTOGRAM_FILE=metrics/unfiltered/$1.insert_metrics.pdf \
									   MINIMUM_PCT=0.5")
									   
metrics/unfiltered/$1.oxog_metrics.txt : bam/$1-unfiltered.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(COLLECT_OXOG_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@)")

metrics/unfiltered/$1.probe-A.hs_metrics.txt : bam/$1-unfiltered.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(CALC_HS_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   BAIT_INTERVALS=$$(POOL_A_TARGET_FILE) \
									   TARGET_INTERVALS=$$(POOL_A_TARGET_FILE)")
 												
metrics/unfiltered/$1.probe-B.hs_metrics.txt : bam/$1-unfiltered.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(CALC_HS_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   BAIT_INTERVALS=$$(POOL_B_TARGET_FILE) \
									   TARGET_INTERVALS=$$(POOL_B_TARGET_FILE)")

endef
$(foreach sample,$(SAMPLES),\
 		$(eval $(call picard-metrics-unfiltered,$(sample))))

define picard-metrics-simplex
metrics/simplex/$1.idx_stats.txt : bam/$1-simplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(BAM_INDEX) \
									   INPUT=$$(<) \
									   > $$(@)")
									   
metrics/simplex/$1.aln_metrics.txt : bam/$1-simplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(COLLECT_ALIGNMENT_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@)")
									   
metrics/simplex/$1.insert_metrics.txt : bam/$1-simplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(COLLECT_INSERT_METRICS) \
 									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   HISTOGRAM_FILE=metrics/simplex/$1.insert_metrics.pdf \
									   MINIMUM_PCT=0.5")
									   
metrics/simplex/$1.oxog_metrics.txt : bam/$1-simplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(COLLECT_OXOG_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@)")

metrics/simplex/$1.probe-A.hs_metrics.txt : bam/$1-simplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(CALC_HS_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   BAIT_INTERVALS=$$(POOL_A_TARGET_FILE) \
									   TARGET_INTERVALS=$$(POOL_A_TARGET_FILE)")
 												
metrics/simplex/$1.probe-B.hs_metrics.txt : bam/$1-simplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(CALC_HS_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   BAIT_INTERVALS=$$(POOL_B_TARGET_FILE) \
									   TARGET_INTERVALS=$$(POOL_B_TARGET_FILE)")

endef
$(foreach sample,$(SAMPLES),\
 		$(eval $(call picard-metrics-simplex,$(sample))))

define picard-metrics-duplex
metrics/duplex/$1.idx_stats.txt : bam/$1-duplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(BAM_INDEX) \
									   INPUT=$$(<) \
									   > $$(@)")
									   
metrics/duplex/$1.aln_metrics.txt : bam/$1-duplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(COLLECT_ALIGNMENT_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@)")
									   
metrics/duplex/$1.insert_metrics.txt : bam/$1-duplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(COLLECT_INSERT_METRICS) \
 									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   HISTOGRAM_FILE=metrics/duplex/$1.insert_metrics.pdf \
									   MINIMUM_PCT=0.5")
									   
metrics/duplex/$1.oxog_metrics.txt : bam/$1-duplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(COLLECT_OXOG_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@)")

metrics/duplex/$1.probe-A.hs_metrics.txt : bam/$1-duplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(CALC_HS_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   BAIT_INTERVALS=$$(POOL_A_TARGET_FILE) \
									   TARGET_INTERVALS=$$(POOL_A_TARGET_FILE)")
 												
metrics/duplex/$1.probe-B.hs_metrics.txt : bam/$1-duplex.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   $$(CALC_HS_METRICS) \
									   REFERENCE_SEQUENCE=$$(REF_FASTA) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   BAIT_INTERVALS=$$(POOL_B_TARGET_FILE) \
									   TARGET_INTERVALS=$$(POOL_B_TARGET_FILE)")

endef
$(foreach sample,$(SAMPLES),\
 		$(eval $(call picard-metrics-duplex,$(sample))))

metrics/summary/umi_frequencies.tsv : $(wildcard marianas/$(SAMPLES)/umi-frequencies.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/umi_metrics.R --type 1 --sample_names '$(SAMPLES)'")
	
metrics/summary/umi_composite.tsv : $(wildcard marianas/$(SAMPLES)/composite-umi-frequencies.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/umi_metrics.R --type 2 --sample_names '$(SAMPLES)'")

metrics/summary/umi_families.tsv : $(wildcard marianas/$(SAMPLES)/family-sizes.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/umi_metrics.R --type 3 --sample_names '$(SAMPLES)'")

metrics/standard/metrics_idx.tsv : $(wildcard metrics/standard/$(SAMPLES).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 1 --sample_names '$(SAMPLES)'")
									  
metrics/unfiltered/metrics_idx.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 7 --sample_names '$(SAMPLES)'")

metrics/duplex/metrics_idx.tsv : $(wildcard metrics/duplex/$(SAMPLES).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 12 --sample_names '$(SAMPLES)'")

metrics/simplex/metrics_idx.tsv : $(wildcard metrics/simplex/$(SAMPLES).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 17 --sample_names '$(SAMPLES)'")
									  
metrics/standard/metrics_aln.tsv : $(wildcard metrics/standard/$(SAMPLES).aln_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 2 --sample_names '$(SAMPLES)'")
									  
metrics/unfiltered/metrics_aln.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).aln_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 8 --sample_names '$(SAMPLES)'")
									  
metrics/duplex/metrics_aln.tsv : $(wildcard metrics/duplex/$(SAMPLES).aln_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 13 --sample_names '$(SAMPLES)'")
									  
metrics/simplex/metrics_aln.tsv : $(wildcard metrics/simplex/$(SAMPLES).aln_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 18 --sample_names '$(SAMPLES)'")
									  
metrics/standard/metrics_insert.tsv : $(wildcard metrics/standard/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 3 --sample_names '$(SAMPLES)'")

metrics/unfiltered/metrics_insert.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 9 --sample_names '$(SAMPLES)'")

metrics/duplex/metrics_insert.tsv : $(wildcard metrics/duplex/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 14 --sample_names '$(SAMPLES)'")

metrics/simplex/metrics_insert.tsv : $(wildcard metrics/simplex/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 19 --sample_names '$(SAMPLES)'")

metrics/standard/metrics_insert_distribution.tsv : $(wildcard metrics/standard/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 4 --sample_names '$(SAMPLES)'")
									   
metrics/unfiltered/metrics_insert_distribution.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 10 --sample_names '$(SAMPLES)'")

metrics/duplex/metrics_insert_distribution.tsv : $(wildcard metrics/duplex/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 15 --sample_names '$(SAMPLES)'")

metrics/simplex/metrics_insert_distribution.tsv : $(wildcard metrics/simplex/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 20 --sample_names '$(SAMPLES)'")
									   
metrics/standard/metrics_hs.tsv : $(wildcard metrics/standard/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/standard/$(SAMPLES).probe-B.hs_metrics.txt) $(wildcard metrics/standard/$(SAMPLES).probe-A.hs_metrics-nodedup.txt) $(wildcard metrics/standard/$(SAMPLES).probe-B.hs_metrics-nodedup.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 6 --sample_names '$(SAMPLES)'")
	
metrics/unfiltered/metrics_hs.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/unfiltered/$(SAMPLES).probe-B.hs_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 11 --sample_names '$(SAMPLES)'")
	
metrics/duplex/metrics_hs.tsv : $(wildcard metrics/duplex/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/duplex/$(SAMPLES).probe-B.hs_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 16 --sample_names '$(SAMPLES)'")
	
metrics/simplex/metrics_hs.tsv : $(wildcard metrics/simplex/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/simplex/$(SAMPLES).probe-B.hs_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 21 --sample_names '$(SAMPLES)'")
									  
metrics/standard/metrics_oxog.tsv : $(wildcard metrics/standard/$(SAMPLES).oxog_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 5 --sample_names '$(SAMPLES)'")
									  
metrics/summary/metrics_idx.tsv : metrics/standard/metrics_idx.tsv metrics/unfiltered/metrics_idx.tsv metrics/duplex/metrics_idx.tsv metrics/simplex/metrics_idx.tsv
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 22")
		
metrics/summary/metrics_aln.tsv : metrics/standard/metrics_aln.tsv metrics/unfiltered/metrics_aln.tsv metrics/duplex/metrics_aln.tsv metrics/simplex/metrics_aln.tsv
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 23")
	
metrics/summary/metrics_insert.tsv : metrics/standard/metrics_insert.tsv metrics/unfiltered/metrics_insert.tsv metrics/duplex/metrics_insert.tsv metrics/simplex/metrics_insert.tsv
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 24")
	
metrics/summary/metrics_insert_distribution.tsv : metrics/standard/metrics_insert_distribution.tsv metrics/unfiltered/metrics_insert_distribution.tsv metrics/duplex/metrics_insert_distribution.tsv metrics/simplex/metrics_insert_distribution.tsv
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 25")
	
metrics/summary/metrics_hs.tsv : metrics/standard/metrics_hs.tsv metrics/unfiltered/metrics_hs.tsv metrics/duplex/metrics_hs.tsv metrics/simplex/metrics_hs.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 26")

metrics/summary/metrics_ts.tsv : $(wildcard metrics/standard/$(SAMPLES).A.ontarget.txt) $(wildcard metrics/standard/$(SAMPLES).B.ontarget.txt) $(wildcard metrics/standard/$(SAMPLES).AB.offtarget.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/interval_metrics.R --metric_type 27 --sample_names '$(SAMPLES)'")									  

..DUMMY := $(shell mkdir -p version; \
			 $(BWA) &> version/tmp.txt; \
			 head -3 version/tmp.txt | tail -2 > version/msk_access.txt; \
			 rm version/tmp.txt; \
			 $(SAMTOOLS) --version >> version/msk_access.txt; \
			 echo "gatk3" >> version/msk_access.txt; \
			 $(GATK) --version >> version/msk_access.txt; \
			 echo "picard" >> version/msk_access.txt; \
			 $(PICARD) MarkDuplicates --version &>> version/msk_access.txt; \
			 R --version >> version.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: msk_access
