include modules/Makefile.inc

LOGDIR ?= log/fusion_catcher.$(NOW)
PHONY += fusion_catcher

FUSION_CATCHER_EXE = $(HOME)/share/usr/src/fusioncatcher/bin/fusioncatcher
FUSION_CATCHER_OPTS = -p 8 -d $(HOME)/share/usr/src/fusioncatcher/data/human_v90

fusion_catcher : $(foreach sample,$(SAMPLES),fusion_catcher/$(sample)/$(sample).1.fastq.gz) \
		 		 $(foreach sample,$(SAMPLES),fusion_catcher/$(sample)/$(sample).2.fastq.gz) \
		 		 $(foreach sample,$(SAMPLES),fusion_catcher/$(sample)/out/taskcomplete)

define fusion-catcher
fusion_catcher/$1/$1.1.fastq.gz : fastq/$1.1.fastq.gz
	$$(call RUN,-c -s 2G -m 4G,"set -o pipefail && \
								mkdir -p fusion_catcher/$1 && \
								cp fastq/$1.1.fastq.gz fusion_catcher/$1/$1.1.fastq.gz")
								
fusion_catcher/$1/$1.2.fastq.gz : fastq/$1.2.fastq.gz
	$$(call RUN,-c -s 2G -m 4G,"set -o pipefail && \
								mkdir -p fusion_catcher/$1 && \
								cp fastq/$1.2.fastq.gz fusion_catcher/$1/$1.2.fastq.gz")
								
fusion_catcher/$1/out/taskcomplete : fusion_catcher/$1/$1.1.fastq.gz fusion_catcher/$1/$1.2.fastq.gz
	$$(call RUN,-c -n 8 -s 2G -m 4G,"set -o pipefail && \
									 mkdir -p fusion_catcher/$1/out && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/bin:$(PATH) && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/tools/bowtie-1.1.2-linux-x86_64:$(PATH) && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/tools/bowtie2-2.2.9-linux-x86_64:$(PATH) && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/tools/bwa/:$(PATH) && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/tools/blat:$(PATH) && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/tools/STAR-2.5.2b/source/:$(PATH) && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/tools/liftover:$(PATH) && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/tools/1.2-r101c:$(PATH) && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/tools/sratoolkit/bin:$(PATH) && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/tools/velvet/:$(PATH) && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/tools/fatotwobit/:$(PATH) && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/tools/lzop/src/:$(PATH) && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/tools/coreutils/src/:$(PATH) && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/tools/pigz/:$(PATH) && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/tools/samtools/:$(PATH) && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/tools/BBMap_37.28/:$(PATH) && \
									 export PATH=$(HOME)/share/usr/src/fusioncatcher/tools/picard/:$(PATH) && \
									 $(FUSION_CATCHER_EXE) $(FUSION_CATCHER_OPTS) -i fusion_catcher/$1 -o fusion_catcher/$1/out && \
									 echo $1 > fusion_catcher/$1/out/taskcomplete")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fusion-catcher,$(sample))))
					

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
