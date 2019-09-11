ifneq ("$(wildcard config.inc)", "")
	include config.inc
endif
ifneq ("$(wildcard project_config.inc)", "")
	include project_config.inc
endif
include modules/config.inc

export

NUM_ATTEMPTS ?= 10
NOW := $(shell date +"%F")
MAKELOG = log/$(@).$(NOW).log

USE_CLUSTER ?= true
QMAKE = modules/scripts/qmake.pl -n $@.$(NOW) $(if $(SLACK_CHANNEL),-c $(SLACK_CHANNEL)) -r $(NUM_ATTEMPTS) -m -s -- make
NUM_JOBS ?= 100

define RUN_QMAKE
$(QMAKE) -e -f $1 -j $2 $(TARGET) && \
	mkdir -p completed_tasks && \
	touch completed_tasks/$@
endef

RUN_MAKE = $(if $(findstring false,$(USE_CLUSTER))$(findstring n,$(MAKEFLAGS)),+$(MAKE) -f $1,$(call RUN_QMAKE,$1,$(NUM_JOBS)))


#==================================================
# Aligners
#==================================================

TARGETS += bwamem
bwamem :
	$(call RUN_MAKE,modules/aligners/bwamem_aligner.mk)

TARGETS += bwa
bwa : NUM_ATTEMPTS = 50
bwa :
	$(call RUN_MAKE,modules/aligners/bwa_aligner.mk)

TARGETS += bowtie
bowtie : NUM_ATTEMPTS = 50
bowtie :
	$(call RUN_MAKE,modules/aligners/bowtie_aligner.mk)

TARGETS += tmap
tmap : NUM_ATTEMPTS = 50
tmap :
	$(call RUN_MAKE,modules/aligners/tmap_aligner.mk)

TARGETS += hisat
hisat : 
	$(call RUN_MAKE,modules/aligners/hisat_aligner.mk)

TARGETS += tophat
tophat : 
	$(call RUN_MAKE,modules/aligners/tophat_aligner.mk)

TARGETS += star
star:
	$(call RUN_MAKE,modules/aligners/star_aligner.mk)

TARGETS += starfusion
starfusion:
	$(call RUN_MAKE,modules/aligners/starfusion_aligner.mk)
	
TARGETS += blast
blast :
	$(call RUN_MAKE,modules/aligners/blast_aligner.mk)
	
TARGETS += msk_access
msk_access :
	$(call RUN_MAKE,modules/test/workflows/msk_access.mk)


#==================================================
# BAM file processing
#==================================================

TARGETS += fix_bam
fix_bam :
	$(call RUN_MAKE,modules/bam_tools/fix_bam.mk)

TARGETS += fix_rg
fix_rg :
	$(call RUN_MAKE,modules/bam_tools/fix_rg.mk)
	
TARGETS += process_bam
process_bam : 
	$(call RUN_MAKE,modules/bam_tools/process_bam.mk)

TARGETS += merge_bam
merge_bam :
	$(call RUN_MAKE,modules/bam_tools/merge_bam.mk)
	
TARGETS += extract_unmapped
extract_unmapped :
	$(call RUN_MAKE,modules/bam_tools/extract_unmapped.mk)
	
TARGETS += bam_to_fasta
bam_to_fasta :
	$(call RUN_MAKE,modules/bam_tools/bam_to_fasta.mk)

	
#==================================================
# FASTQ file processing
#==================================================

TARGETS += extract_fastq
extract_fastq :
	$(call RUN_MAKE,modules/fastq_tools/extract_fastq.mk)

TARGETS += merge_fastq
merge_fastq : 
	$(call RUN_MAKE,modules/fastq_tools/merge_fastq.mk)

TARGETS += merge_split_fastq
merge_split_fastq :
	$(call RUN_MAKE,modules/fastq_tools/merge_split_fastq.mk)
	
	
#==================================================
# Somatic variant callers
#==================================================

TARGETS += mutect
mutect :
	$(call RUN_MAKE,modules/variant_callers/somatic/mutect.mk)
	
TARGETS += varscan
varscan :
	$(call RUN_MAKE,modules/variant_callers/somatic/varscan.mk)
	
TARGETS += strelka
strelka :
	$(call RUN_MAKE,modules/variant_callers/somatic/strelka.mk)
	
TARGETS += scalpel
scalpel :
	$(call RUN_MAKE,modules/variant_callers/somatic/scalpel.mk)
    
TARGETS += lancet
lancet :
	$(call RUN_MAKE,modules/variant_callers/somatic/lancet.mk)
	
TARGETS += platypus
platypus:
	$(call RUN_MAKE,modules/variant_callers/somatic/platypus.mk)
	
	
TARGETS += somatic_sniper
somatic_sniper :
	$(call RUN_MAKE,modules/variant_callers/somatic/somatic_sniper.mk)

TARGETS += tvc_tn
tvc_tn:
	$(call RUN_MAKE,modules/variant_callers/somatic/tvc_tn.mk)
	
	
#==================================================
# Specialty somatic variant callers
#==================================================
	
TARGETS += hla_polysolver
hla_polysolver :
	$(call RUN_MAKE,modules/variant_callers/somatic/hla_polysolver.mk)
	
TARGETS += msi_sensor
msi_sensor :
	$(call RUN_MAKE,modules/variant_callers/somatic/msi_sensor.mk)
	

#==================================================
# Single sample variant callers
#==================================================

TARGETS += tvc
tvc:
	$(call RUN_MAKE,modules/variant_callers/tvc.mk)

TARGETS += gatk
gatk : 
	$(call RUN_MAKE,modules/variant_callers/gatk.mk)
	
TARGETS += haplotype_caller
haplotype_caller : 
	$(call RUN_MAKE,modules/variant_callers/haplotype_caller.mk)
	
TARGETS += samtools_het
samtools_het :
	$(call RUN_MAKE,modules/variant_callers/samtoolsHet.mk)

TARGETS += hotspot
hotspot: 
	$(call RUN_MAKE,modules/variant_callers/hotspot.mk)
	
TARGETS += genotype_hotspot
genotype_hotspot:
	$(call RUN_MAKE,modules/variant_callers/genotypehotspots.mk)
	
TARGETS += genotype_pdx
genotype_pdx:
	$(call RUN_MAKE,modules/variant_callers/genotypepdx.mk)
	
TARGETS += sufam
sufam:
	$(call RUN_MAKE,modules/variant_callers/sufamsampleset.mk)
	
TARGETS += sufam_summary
sufam_summary:
	$(call RUN_MAKE,modules/variant_callers/sufammultisample.mk)


#==================================================
# Copy number aberration callers
#==================================================

TARGETS += facets
facets :
	$(call RUN_MAKE,modules/copy_number/facets.mk)
	
TARGETS += ascat
ascat :
	$(call RUN_MAKE,modules/copy_number/ascat.mk)
	
TARGETS += cnvkit
cnvkit :
	$(call RUN_MAKE,modules/copy_number/cnvkit.mk)

TARGETS += titan
titan :
	$(call RUN_MAKE,modules/copy_number/titan.mk)

TARGETS += varscan_cnv
varscan_cnv :
	$(call RUN_MAKE,modules/copy_number/varscan.mk)

TARGETS += hmmcopy
hmmcopy :
	$(call RUN_MAKE,modules/copy_number/hmmcopy.mk)

TARGETS += gistic
gistic :
	$(call RUN_MAKE,modules/copy_number/gistic.mk)
	
TARGETS += snp6
snp6 :
	$(call RUN_MAKE,modules/copy_number/snp6.mk)
	
TARGETS += qdnaseq
qdnaseq :
	$(call RUN_MAKE,modules/test/workflows/qdna_seq.mk)
	
TARGETS += copynumber_summary
copynumber_summary:
	$(call RUN_MAKE,modules/test/workflows/copynumber_summary.mk)
	

#==================================================
# Structural variant callers
#==================================================

TARGETS += star_fusion
star_fusion:
	$(call RUN_MAKE,modules/sv_callers/starFusion.mk)

TARGETS += tophat_fusion
tophat_fusion : 
	$(call RUN_MAKE,modules/sv_callers/tophatFusion.mk)

TARGETS += manta_rnaseq
manta_rnaseq :
	$(call RUN_MAKE,modules/sv_callers/mantaRnaseq.mk)

TARGETS += manta
manta :
	$(call RUN_MAKE,modules/sv_callers/manta.mk)

TARGETS += mantaTN
mantaTN :
	$(call RUN_MAKE,modules/sv_callers/mantaTN.mk)

TARGETS += brass
brass :
	$(call RUN_MAKE,modules/sv_callers/brass.mk)

TARGETS += integrate_rnaseq
integrate_rnaseq :
	$(call RUN_MAKE,modules/sv_callers/integrateRnaseq.mk)

TARGETS += integrate
integrate :
	$(call RUN_MAKE,modules/sv_callers/integrate.mk)

TARGETS += defuse
defuse :
	$(call RUN_MAKE,modules/sv_callers/defuse.mk)

NUM_CHIMSCAN_JOBS ?= 5
TARGETS += chimscan
chimscan :
	$(call RUN_MAKE_J,modules/sv_callers/chimerascan.mk,$(NUM_CHIMSCAN_JOBS))

TARGETS += oncofuse
oncofuse :
	$(call RUN_MAKE,modules/sv_callers/oncofuse.mk)

TARGETS += lumpy
lumpy :
	$(call RUN_MAKE,modules/sv_callers/lumpy.mk)

TARGETS += hydra
hydra :
	$(call RUN_MAKE,modules/sv_callers/hydra.mk)

TARGETS += nfuse_wgss_wtss
nfuse_wgss_wtss :
	$(call RUN_MAKE,modules/sv_callers/nfuseWGSSWTSS.mk)

TARGETS += soapfuse
soapfuse :
	$(call RUN_MAKE,modules/sv_callers/soapFuse.mk)

TARGETS += mapsplice
mapsplice :
	$(call RUN_MAKE,modules/sv_callers/mapsplice.mk)

TARGETS += fusioncatcher
fusioncatcher :
	$(call RUN_MAKE,modules/sv_callers/fusioncatcher.mk)

TARGETS += crest
crest :
	$(call RUN_MAKE,modules/sv_callers/crest.mk)

TARGETS += delly
delly :
	$(call RUN_MAKE,modules/sv_callers/delly.mk)


#==================================================
# Quality control
#==================================================

TARGETS += bam_metrics
bam_metrics :
	$(call RUN_MAKE,modules/qc/bam_metrics.mk)

TARGETS += bam_interval_metrics
bam_interval_metrics :
	$(call RUN_MAKE,modules/qc/bam_interval_metrics.mk)

TARGETS += rnaseq_metrics
rnaseq_metrics :
	$(call RUN_MAKE,modules/qc/rnaseq_metrics.mk)

TARGETS += fastqc
fastqc :
	$(call RUN_MAKE,modules/qc/fastqc.mk)

TARGETS += interval_qc
interval_qc :
	$(call RUN_MAKE,modules/qc/interval_bam_qc.mk)

TARGETS += rseqc
rseqc :
	$(call RUN_MAKE,modules/qc/rseqc.mk)

TARGETS += qualimap
qualimap :
	$(call RUN_MAKE,modules/qc/qualimap_qc.mk)

TARGETS += bam_stats
bam_stats :
	$(call RUN_MAKE,modules/qc/bam_stats.mk)
	
TARGETS += cluster_samples
cluster_samples :
	$(call RUN_MAKE,modules/qc/cluster_samples.mk)


#==================================================
# RNA sequencing
#==================================================

TARGETS += cufflinks
cufflinks : 
	$(call RUN_MAKE,modules/rnaseq/cufflinks.mk)

TARGETS += sum_reads
sum_reads :
	$(call RUN_MAKE,modules/rnaseq/sumRNASeqReads.mk)

TARGETS += exon_counts
exon_counts :
	$(call RUN_MAKE,modules/rnaseq/dexseq.mk)
	

#==================================================
# Ploidy
#==================================================

TARGETS += pyloh
pyloh :
	$(call RUN_MAKE,modules/ploidy/pyloh.mk)


#==================================================
# Clonality
#==================================================

TARGETS += clonehd
clonehd :
	$(call RUN_MAKE,modules/clonality/clonehd.mk)

TARGETS += absolute
absolute :
	$(call RUN_MAKE,modules/clonality/absolute.mk)
	
TARGETS += sspyclone
sspyclone :
	$(call RUN_MAKE,modules/clonality/sspyclone.mk)

TARGETS += mspyclone
mspyclone :
	$(call RUN_MAKE,modules/clonality/mspyclone.mk)
	

#==================================================
# Mutational signatures
#==================================================

TARGETS += emu
emu :
	$(call RUN_MAKE,modules/signatures/emu.mk)
	
TARGETS += mut_sig
mut_sig :
	$(call RUN_MAKE,modules/signatures/mut_sig.mk)
	
TARGETS += deconstruct_sigs
deconstruct_sigs :
	$(call RUN_MAKE,modules/signatures/deconstruct_sigs.mk)


#==================================================
# Miscellaneous
#==================================================

TARGETS += contest
contest :
	$(call RUN_MAKE,modules/contamination/contest.mk)

TARGETS += virus_detection_bowtie2
virus_detection_bowtie2 :
	$(call RUN_MAKE,modules/virus/virus_detection_bowtie2.mk)
	
TARGETS += viral_detection
viral_detection:
	$(call RUN_MAKE,modules/test/workflows/viral_detection.mk)
	
TARGETS += krona_classify
krona_classify :
	$(call RUN_MAKE,modules/virus/krona_classify.mk)
	
TARGETS += fetch_impact
fetch_impact :
	$(call RUN_MAKE,modules/test/workflows/fetch_impact.mk)


#==================================================
# Phylogeny
#==================================================

TARGETS += medicc
medicc :
	$(call RUN_MAKE,modules/test/workflows/medicc.mk)
	
TARGETS += pratchet
pratchet :
	$(call RUN_MAKE,modules/test/workflows/pratchet.mk)


#==================================================
# Reports
#==================================================

TARGETS += recurrent_mutations
recurrent_mutations :
	$(call RUN_MAKE,modules/recurrent_mutations/report.mk)
	
TARGETS += genome_summary
genome_summary :
	$(call RUN_MAKE,modules/summary/genomesummary.mk)

TARGETS += mutation_summary
mutation_summary :
	$(call RUN_MAKE,modules/summary/mutationsummary.mk)
	
TARGETS += cravat_summary
cravat_summary :
	$(call RUN_MAKE,modules/summary/cravat_summary.mk)


#==================================================
# Annotations
#==================================================

TARGETS += ann_ext_vcf
ann_ext_vcf: 
	$(call RUN_MAKE,modules/vcf_tools/annotateExtVcf.mk)

TARGETS += ann_somatic_vcf
ann_somatic_vcf: 
	$(call RUN_MAKE,modules/vcf_tools/annotateSomaticVcf.mk)

TARGETS += ann_vcf
ann_vcf: 
	$(call RUN_MAKE,modules/vcf_tools/annotateVcf.mk)
	
TARGETS += cravat_annotation
cravat_annotation :
	$(call RUN_MAKE,modules/test/workflows/cravat_annotation.mk)
	
TARGETS += cravat_annotate
cravat_annotate :
	$(call RUN_MAKE,modules/vcf_tools/cravat_annotation.mk)


#==================================================
# Alpha testing
#==================================================

TARGETS += hotspot_summary
hotspot_summary:
	$(MAKE) -f modules/variant_callers/genotypehotspots.mk -j $(NUM_JOBS)
	$(call RUN_MAKE,modules/summary/hotspotsummary.mk)
	
TARGETS += cnvkit_coverage
cnvkit_coverage :
	$(call RUN_MAKE,modules/copy_number/cnvkitcoverage.mk)
	
TARGETS += cnvkit_reference
cnvkit_reference :
	$(call RUN_MAKE,modules/copy_number/cnvkitreference.mk)
	
TARGETS += cnvkit_fix
cnvkit_fix :
	$(call RUN_MAKE,modules/copy_number/cnvkitfix.mk)

TARGETS += cnvkit_plot
cnvkit_plot :
	$(call RUN_MAKE,modules/copy_number/cnvkitplot.mk)
	
TARGETS += cnvkit_heatmap
cnvkit_heatmap :
	$(call RUN_MAKE,modules/copy_number/cnvkitheatmap.mk)
	
TARGETS += cnvkit_pca
cnvkit_pca :
	$(call RUN_MAKE,modules/copy_number/cnvkitprcomp.mk)
	
TARGETS += cnvkit_qc
cnvkit_qc :
	$(call RUN_MAKE,modules/copy_number/cnvkitqc.mk)
	

#==================================================
# Beta testing
#==================================================

TARGETS += somatic_indels
somatic_indels:
	$(call RUN_MAKE,modules/test/workflows/somatic_indels.mk)
	
TARGETS += somatic_variants
somatic_variants:
	$(call RUN_MAKE,modules/test/workflows/somatic_variants.mk)
	
TARGETS += fgbio_access
fgbio_access :
	$(call RUN_MAKE,modules/test/workflows/fgbio_access.mk)
	
TARGETS += cnv_access
cnv_access :
	$(call RUN_MAKE,modules/test/workflows/cnv_access.mk)
	

.PHONY : $(TARGETS)
