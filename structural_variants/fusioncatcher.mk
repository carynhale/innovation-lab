include modules/Makefile.inc

LOGDIR ?= log/fusion_catcher.$(NOW)
PHONY += fusion_catcher

PATH=$(HOME)/share/usr/src/fusioncatcher/bin:$PATH
PATH=$(HOME)/share/usr/src/fusioncatcher/tools/bowtie-1.1.2-linux-x86_64:$PATH
PATH=$(HOME)/share/usr/src/fusioncatcher/tools/bowtie2-2.2.9-linux-x86_64:$PATH
PATH=$(HOME)/share/usr/src/fusioncatcher/tools/bwa/:$PATH
PATH=$(HOME)/share/usr/src/fusioncatcher/tools/blat:$PATH
PATH=$(HOME)/share/usr/src/fusioncatcher/tools/STAR-2.5.2b/source/:$PATH
PATH=$(HOME)/share/usr/src/fusioncatcher/tools/liftover:$PATH
PATH=$(HOME)/share/usr/src/fusioncatcher/tools/1.2-r101c:$PATH
PATH=$(HOME)/share/usr/src/fusioncatcher/tools/sratoolkit/bin:$PATH
PATH=$(HOME)/share/usr/src/fusioncatcher/tools/velvet/:$PATH
PATH=$(HOME)/share/usr/src/fusioncatcher/tools/fatotwobit/:$PATH
PATH=$(HOME)/share/usr/src/fusioncatcher/tools/lzop/src/:$PATH
PATH=$(HOME)/share/usr/src/fusioncatcher/tools/coreutils/src/:$PATH
PATH=$(HOME)/share/usr/src/fusioncatcher/tools/pigz/:$PATH
PATH=$(HOME)/share/usr/src/fusioncatcher/tools/samtools/:$PATH
PATH=$(HOME)/share/usr/src/fusioncatcher/tools/BBMap_37.28/:$PATH
PATH=$(HOME)/share/usr/src/fusioncatcher/tools/picard/:$PATH
FUSION_CATCHER_EXE = $(HOME)/share/usr/fusioncatcher/bin/fusioncatcher
FUSIONCATCHER_OPTS = -d $(HOME)/share/usr/fusioncatcher/data/current --extract-buffer-size=35000000000

#fusioncatcher -d ~/share/usr/src/fusioncatcher/data/human_v90 -i rawdata/PITT_0392/Sample_MCM101T_IGO_04835_J_1/ -o
#fusion_catcher/%/%.taskcomplete : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
#	$(call RUN,-n 8 -s 1G -m 4G,"$(FUSIONCATCHER) $(FUSIONCATCHER_OPTS) -p 8 -o $(@D)/$* -i $<$(,)$(<<) && \
#								 touch fusion-ctahcer")
#
#	

fusion_catcher : $(foreach sample,$(SAMPLES),fusion_catcher/$(sample)/$(sample).1.fastq.gz) \
		 		 $(foreach sample,$(SAMPLES),fusion_catcher/$(sample)/$(sample).2.fastq.gz)

define fusion-catcher
fusion_catcher/$1/$1.1.fastq.gz : fastq/$1.1.fastq.gz
	$$(call RUN,-c -s 2G -m 4G,"mkdir -p fusion_catcher/$$(*) && \
								cp fastq/$$(*).1.fastq.gz fusion_catcher/$$(*)/$$(*).1.fastq.gz")
								
fusion_catcher/$1/$1.2.fastq.gz : fastq/$1.2.fastq.gz
	$$(call RUN,-c -s 2G -m 4G,"mkdir -p fusion_catcher/$$(*) && \
								cp fastq/$$(*).2.fastq.gz fusion_catcher/$$(*)/$$(*).2.fastq.gz")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fusion-catcher,$(sample))))
					

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
