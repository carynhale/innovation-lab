ifndef GATK_INC
GATK_MEM = 18G
GATK_CMD = $(call GATK_OPT,$(GATK_MEM))
GATK_OPT = $(GATK) -Xmx$(1) -S LENIENT

endif
GATK_INC = true
