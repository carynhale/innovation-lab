ifneq ("$(wildcard config.inc)", "")
	include config.inc
endif
ifneq ("$(wildcard project_config.inc)", "")
	include project_config.inc
endif
include innovation-lab/config/config.inc

export

NUM_ATTEMPTS ?= 10
NOW := $(shell date +"%F")
MAKELOG = log/$(@).$(NOW).log

USE_CLUSTER ?= true
QMAKE = innovation-lab/dodo-cloning-kit/runtime/qmake.pl -n $@.$(NOW) $(if $(SLACK_CHANNEL),-c $(SLACK_CHANNEL)) -r $(NUM_ATTEMPTS) -m -s -- make
NUM_JOBS ?= 50

define RUN_QMAKE
$(QMAKE) -e -f $1 -j $2 $(TARGET) && \
	mkdir -p completed_tasks && \
	touch completed_tasks/$@
endef

RUN_MAKE = $(if $(findstring false,$(USE_CLUSTER))$(findstring n,$(MAKEFLAGS)),+$(MAKE) -f $1,$(call RUN_QMAKE,$1,$(NUM_JOBS)))

#==================================================
# Workflows
#==================================================

TARGETS += msk_access
msk_access :
	$(call RUN_MAKE,innovation-lab/workflows/msk_access.mk)


#==================================================
# FASTQ / BAM file aligners
#==================================================

TARGETS += bwa_mem
bwa_mem :
	$(call RUN_MAKE,innovation-lab/aligners/bwa_mem.mk)
	

#==================================================
# BAM file processing utilities
#==================================================

TARGETS += fix_bam
fix_bam :
	$(call RUN_MAKE,innovation-lab/bam_tools/fix_bam.mk)

TARGETS += unmapped_fasta
unmapped_fasta :
	$(call RUN_MAKE,innovation-lab/bam_tools/unmapped_fasta.mk)
	
TARGETS += unmapped_bam
unmapped_bam :
	$(call RUN_MAKE,innovation-lab/bam_tools/unmapped_bam.mk)

TARGETS += merge_bam
merge_bam :
	$(call RUN_MAKE,innovation-lab/bam_tools/merge_bam.mk)
	
	
#==================================================
# FASTQ file processing utilities
#==================================================

TARGETS += extract_fastq
extract_fastq :
	$(call RUN_MAKE,innovation-lab/fastq_tools/extract_fastq.mk)

TARGETS += merge_fastq
merge_fastq : 
	$(call RUN_MAKE,innovation-lab/fastq_tools/merge_fastq.mk)

TARGETS += merge_split_fastq
merge_split_fastq :
	$(call RUN_MAKE,innovation-lab/fastq_tools/merge_split_fastq.mk)
	
	
#==================================================
# Somatic variant callers
#==================================================

TARGETS += mutect
mutect :
	$(call RUN_MAKE,innovation-lab/variant_callers/mutect.mk)
	
TARGETS += varscan
varscan :
	$(call RUN_MAKE,innovation-lab/variant_callers/varscan.mk)
	
TARGETS += strelka
strelka :
	$(call RUN_MAKE,innovation-lab/variant_callers/strelka.mk)
	
TARGETS += scalpel
scalpel :
	$(call RUN_MAKE,innovation-lab/variant_callers/scalpel.mk)
    
TARGETS += lancet
lancet :
	$(call RUN_MAKE,innovation-lab/variant_callers/lancet.mk)
	
TARGETS += platypus
platypus:
	$(call RUN_MAKE,innovation-lab/variant_callers/platypus.mk)
	
	
#==================================================
# Specialty variant callers
#==================================================
	
TARGETS += hla_polysolver
hla_polysolver :
	$(call RUN_MAKE,innovation-lab/variant_callers/somatic/hla_polysolver.mk)
	
TARGETS += msi_sensor
msi_sensor :
	$(call RUN_MAKE,innovation-lab/variant_callers/somatic/msi_sensor.mk)

TARGETS += haplotype_caller
haplotype_caller : 
	$(call RUN_MAKE,innovation-lab/variant_callers/haplotype_caller.mk)
	
TARGETS += genotype_hotspot
genotype_hotspot: 
	$(call RUN_MAKE,innovation-lab/variant_callers/genotype_hotspots.mk)
	
TARGETS += multisample_genotype
multisample_genotype:
	$(call RUN_MAKE,innovation-lab/variant_callers/multisample_genotype.mk)
	
TARGETS += waltz_genotype
waltz_genotype:
	$(call RUN_MAKE,innovation-lab/variant_callers/waltz_genotype.mk)


#==================================================
# Copy number aberration callers
#==================================================

TARGETS += facets
facets :
	$(call RUN_MAKE,innovation-lab/copy_number/facets.mk)
	
TARGETS += ascat
ascat :
	$(call RUN_MAKE,innovation-lab/copy_number/ascat.mk)
	
TARGETS += cnvkit
cnvkit :
	$(call RUN_MAKE,innovation-lab/copy_number/cnvkit.mk)

TARGETS += qdnaseq
qdnaseq :
	$(call RUN_MAKE,innovation-lab/test/workflows/qdna_seq.mk)
	
TARGETS += cnvaccess
cnvaccess :
	$(call RUN_MAKE,innovation-lab/test/workflows/cnv_access.mk)
	

#==================================================
# Structural variant callers
#==================================================

TARGETS += rnaseq_fusion
rnaseq_fusion:
	$(call RUN_MAKE,innovation-lab/test/workflows/rnaseq_fusion.mk)

TARGETS += star_fusion
star_fusion:
	$(call RUN_MAKE,innovation-lab/structural_variants/star_fusion.mk)

TARGETS += tophat_fusion
tophat_fusion : 
	$(call RUN_MAKE,innovation-lab/structural_variants/tophat_fusion.mk)

TARGETS += manta_rnaseq
manta_rnaseq :
	$(call RUN_MAKE,innovation-lab/structural_variants/manta_rnaseq.mk)

TARGETS += manta
manta :
	$(call RUN_MAKE,innovation-lab/structural_variants/manta.mk)

TARGETS += manta_tn
manta_tn :
	$(call RUN_MAKE,innovation-lab/structural_variants/manta_tn.mk)

TARGETS += brass
brass :
	$(call RUN_MAKE,innovation-lab/structural_variants/brass.mk)

TARGETS += integrate
integrate :
	$(call RUN_MAKE,innovation-lab/structural_variants/integrate.mk)
	
TARGETS += integrate_rnaseq
integrate_rnaseq :
	$(call RUN_MAKE,innovation-lab/structural_variants/integrate_rnaseq.mk)

TARGETS += defuse
defuse :
	$(call RUN_MAKE,innovation-lab/structural_variants/defuse.mk)

NUM_CHIMSCAN_JOBS ?= 5
TARGETS += chimscan
chimscan :
	$(call RUN_MAKE_J,innovation-lab/structural_variants/chimerascan.mk,$(NUM_CHIMSCAN_JOBS))

TARGETS += oncofuse
oncofuse :
	$(call RUN_MAKE,innovation-lab/structural_variants/oncofuse.mk)

TARGETS += lumpy
lumpy :
	$(call RUN_MAKE,innovation-lab/structural_variants/lumpy.mk)

TARGETS += hydra
hydra :
	$(call RUN_MAKE,innovation-lab/structural_variants/hydra.mk)

TARGETS += nfuse_wgss_wtss
nfuse_wgss_wtss :
	$(call RUN_MAKE,innovation-lab/structural_variants/nfuseWGSSWTSS.mk)

TARGETS += soap_fuse
soap_fuse :
	$(call RUN_MAKE,innovation-lab/structural_variants/soap_fuse.mk)

TARGETS += map_splice
map_splice :
	$(call RUN_MAKE,innovation-lab/structural_variants/map_splice.mk)

TARGETS += fusion_catcher
fusion_catcher :
	$(call RUN_MAKE,innovation-lab/structural_variants/fusion_catcher.mk)

TARGETS += crest
crest :
	$(call RUN_MAKE,innovation-lab/structural_variants/crest.mk)

TARGETS += delly
delly :
	$(call RUN_MAKE,innovation-lab/structural_variants/delly.mk)


#==================================================
# Quality control
#==================================================

TARGETS += bam_metrics
bam_metrics :
	$(call RUN_MAKE,innovation-lab/qc/bam_metrics.mk)

TARGETS += bam_interval_metrics
bam_interval_metrics :
	$(call RUN_MAKE,innovation-lab/qc/bam_interval_metrics.mk)

TARGETS += rnaseq_metrics
rnaseq_metrics :
	$(call RUN_MAKE,innovation-lab/qc/rnaseq_metrics.mk)

TARGETS += fastqc
fastqc :
	$(call RUN_MAKE,innovation-lab/qc/fastqc.mk)

TARGETS += interval_qc
interval_qc :
	$(call RUN_MAKE,innovation-lab/qc/interval_bam_qc.mk)

TARGETS += rseqc
rseqc :
	$(call RUN_MAKE,innovation-lab/qc/rseqc.mk)

TARGETS += qualimap
qualimap :
	$(call RUN_MAKE,innovation-lab/qc/qualimap_qc.mk)

TARGETS += bam_stats
bam_stats :
	$(call RUN_MAKE,innovation-lab/qc/bam_stats.mk)
	
TARGETS += cluster_samples
cluster_samples :
	$(call RUN_MAKE,innovation-lab/qc/cluster_samples.mk)


#==================================================
# RNA sequencing
#==================================================

TARGETS += cufflinks
cufflinks : 
	$(call RUN_MAKE,innovation-lab/rnaseq/cufflinks.mk)

TARGETS += sum_reads
sum_reads :
	$(call RUN_MAKE,innovation-lab/rnaseq/sumrnaseqreads.mk)

TARGETS += exon_counts
exon_counts :
	$(call RUN_MAKE,innovation-lab/rnaseq/dexseq.mk)
	

#==================================================
# Clonality
#==================================================

TARGETS += absolute
absolute :
	$(call RUN_MAKE,innovation-lab/clonality/absolute.mk)
	
TARGETS += pyclone
pyclone :
	$(call RUN_MAKE,innovation-lab/clonality/pyclone.mk)

#==================================================
# Miscellaneous
#==================================================

TARGETS += viral_detection
viral_detection:
	$(call RUN_MAKE,innovation-lab/test/workflows/viral_detection.mk)
	
TARGETS += krona_classify
krona_classify :
	$(call RUN_MAKE,innovation-lab/virus/krona_classify.mk)
	
TARGETS += fetch_impact
fetch_impact :
	$(call RUN_MAKE,innovation-lab/test/workflows/fetch_impact.mk)


#==================================================
# Reports
#==================================================

TARGETS += recurrent_mutations
recurrent_mutations :
	$(call RUN_MAKE,innovation-lab/recurrent_mutations/report.mk)
	
TARGETS += genome_summary
genome_summary :
	$(call RUN_MAKE,innovation-lab/summary/genomesummary.mk)

TARGETS += mutation_summary
mutation_summary :
	$(call RUN_MAKE,innovation-lab/summary/mutationsummary.mk)
	
TARGETS += cravat_summary
cravat_summary :
	$(call RUN_MAKE,innovation-lab/summary/cravat_summary.mk)
	
TARGETS += copynumber_summary
copynumber_summary:
	$(call RUN_MAKE,innovation-lab/test/workflows/copynumber_summary.mk)
	
TARGETS += somatic_indels
somatic_indels:
	$(call RUN_MAKE,innovation-lab/test/workflows/somatic_indels.mk)
	
TARGETS += somatic_variants
somatic_variants:
	$(call RUN_MAKE,innovation-lab/test/workflows/somatic_variants.mk)
	

#==================================================
# Annotations
#==================================================

TARGETS += ann_ext_vcf
ann_ext_vcf: 
	$(call RUN_MAKE,innovation-lab/vcf_tools/annotateExtVcf.mk)

TARGETS += ann_somatic_vcf
ann_somatic_vcf: 
	$(call RUN_MAKE,innovation-lab/vcf_tools/annotateSomaticVcf.mk)

TARGETS += ann_vcf
ann_vcf: 
	$(call RUN_MAKE,innovation-lab/vcf_tools/annotateVcf.mk)
	
TARGETS += cravat_annotation
cravat_annotation :
	$(call RUN_MAKE,innovation-lab/test/workflows/cravat_annotation.mk)
	
TARGETS += cravat_annotate
cravat_annotate :
	$(call RUN_MAKE,innovation-lab/vcf_tools/cravat_annotation.mk)


#==================================================
# Alpha testing
#==================================================

TARGETS += hotspot_summary
hotspot_summary:
	$(MAKE) -f innovation-lab/variant_callers/genotypehotspots.mk -j $(NUM_JOBS)
	$(call RUN_MAKE,innovation-lab/summary/hotspotsummary.mk)
	
TARGETS += cnvkit_coverage
cnvkit_coverage :
	$(call RUN_MAKE,innovation-lab/copy_number/cnvkitcoverage.mk)
	
TARGETS += cnvkit_reference
cnvkit_reference :
	$(call RUN_MAKE,innovation-lab/copy_number/cnvkitreference.mk)
	
TARGETS += cnvkit_fix
cnvkit_fix :
	$(call RUN_MAKE,innovation-lab/copy_number/cnvkitfix.mk)

TARGETS += cnvkit_plot
cnvkit_plot :
	$(call RUN_MAKE,innovation-lab/copy_number/cnvkitplot.mk)
	
TARGETS += cnvkit_heatmap
cnvkit_heatmap :
	$(call RUN_MAKE,innovation-lab/copy_number/cnvkitheatmap.mk)
	
TARGETS += cnvkit_pca
cnvkit_pca :
	$(call RUN_MAKE,innovation-lab/copy_number/cnvkitprcomp.mk)
	
TARGETS += cnvkit_qc
cnvkit_qc :
	$(call RUN_MAKE,innovation-lab/copy_number/cnvkitqc.mk)
	

TARGETS += fgbio_access
fgbio_access :
	$(call RUN_MAKE,innovation-lab/test/workflows/fgbio_access.mk)
	
.PHONY : $(TARGETS)
