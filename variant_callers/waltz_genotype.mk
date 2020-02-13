include innovation-lab/Makefile.inc
include innovation-lab/config/waltz.inc
include innovation-lab/genome_inc/b37.inc

..DUMMY := $(shell mkdir -p version)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: waltz_genotype

