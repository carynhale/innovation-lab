include modules/Makefile.inc

LOGDIR ?= log/fusion_catcher.$(NOW)
.PHONY += fusion_catcher

export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/bin:$PATH
export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/tools/bowtie-1.1.2-linux-x86_64:$PATH
export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/tools/bowtie2-2.2.9-linux-x86_64:$PATH
export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/tools/bwa/:$PATH
export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/tools/blat:$PATH
export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/tools/STAR-2.5.2b/source/:$PATH
export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/tools/liftover:$PATH
export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/tools/1.2-r101c:$PATH
export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/tools/sratoolkit/bin:$PATH
export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/tools/velvet/:$PATH
export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/tools/fatotwobit/:$PATH
export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/tools/lzop/src/:$PATH
export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/tools/coreutils/src/:$PATH
export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/tools/pigz/:$PATH
export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/tools/samtools/:$PATH
export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/tools/BBMap_37.28/:$PATH
export PATH=/home/dacruzpa/share/usr/src/fusioncatcher/tools/picard/:$PATH

fusion_catcher : $(foreach sample,$(SAMPLES),fusioncatcher/$(sample)/$(sample).taskcomplete)

FUSION_CATCHER_ = $(HOME)/share/usr/fusioncatcher/bin/fusioncatcher
FUSIONCATCHER_OPTS = -d $(HOME)/share/usr/fusioncatcher/data/current --extract-buffer-size=35000000000

define fusion-catcher
fusion_catcher/%/%.taskcomplete : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call RUN,-n 8 -s 1G -m 4G,"$(FUSIONCATCHER) $(FUSIONCATCHER_OPTS) -p 8 -o $(@D)/$* -i $<$(,)$(<<) && \
								 touch fusion-ctahcer")

	
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fusion-catcher,$(sample))))
					

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)


fusioncatcher -d ~/share/usr/src/fusioncatcher/data/human_v90 -i rawdata/PITT_0392/Sample_MCM101T_IGO_04835_J_1/ -o