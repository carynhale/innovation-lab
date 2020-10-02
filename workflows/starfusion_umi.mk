include innovation-lab/Makefile.inc

LOGDIR ?= log/starfusion_umi.$(NOW)

starfusion_umi : $(foreach sample,$(SAMPLES),starfusion/$(sample).dedup/$(sample).1.fastq) \
				 $(foreach sample,$(SAMPLES),starfusion/$(sample).dedup/taskcomplete) \
				 starfusion/summary.dedup.txt
			  
CTAT_LIB ?= $(HOME)/share/lib/ref_files/CTAT_GRCh37/GRCh37_gencode_v19_CTAT_lib_Apr032020/ctat_genome_lib_build_dir/

define starfusion-dedup
starfusion/$1.dedup/$1.1.fastq : bam/$1.bam
	$$(call RUN,-n 4 -s 4G -m 9G,"set -o pipefail && \
								  $$(SAMTOOLS) sort -T bam/$1 -O bam -n -@ 4 -m 6G $$(<) | \
								  $$(SAMTOOLS) fastq -f 1 -1 fusioncatcher/$1.dedup/$1.1.fastq -2 fusioncatcher/$1.dedup/$1.2.fastq")
								  
starfusion/$1.dedup/taskcomplete : starfusion/$1.dedup/$1.1.fastq
	$$(call RUN,-c -n 20 -s 1G -m 2G -v $(STARFUSION_ENV) -w 36:00:00,"set -o pipefail && \
																	   $$(STAR_FUSION) \
																	   --left_fq starfusion/$1.dedup/$1.1.fastq \
																	   --right_fq starfusion/$1.dedup/$1.2.fastq \
																	   --CPU 20 \
																	   --output_dir starfusion/$1.dedup \
																	   --genome_lib_dir $$(CTAT_LIB) && \
																	   echo $1 > starfusion/$1.dedup/taskcomplete")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call starfusion-dedup,$(sample))))
		
starfusion/summary.dedup.txt : $(foreach sample,$(SAMPLES),starfusion/$(sample).dedup/taskcomplete)
	echo "FusionName	JunctionReadCount	SpanningFragCount	SpliceType	LeftGene	LeftBreakpoint	RightGene	RightBreakpoint	LargeAnchorSupport	FFPM	LeftBreakDinuc	LeftBreakEntropy	RightBreakDinuc	RightBreakEntropy	annots	SampleName" > starfusion/summary.dedup.txt; \
	for i in $(SAMPLES); do \
		sed -e "1d" starfusion/$$i.dedup/star-fusion.fusion_predictions.abridged.tsv | sed "s/$$/\t$$i/" >> starfusion/summary.dedup.txt; \
	done
	

..DUMMY := $(shell mkdir -p version; \
			 ~/share/usr/env/starfusion-1.6.0/bin/STAR-Fusion --help &> version/starfusion_umi.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: starfusion_umi
