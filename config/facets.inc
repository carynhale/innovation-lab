FACETS_HET_THRESHOLD ?= 0.25
FACETS_MIN_NHET ?= 25
FACETS_SNP_NBHD ?= 250
FACETS_NDEPTH_MAX ?= 1000
FACETS_PRE_CVAL ?= 50
FACETS_CVAL1 ?= 150
FACETS_CVAL2 ?= 50
FACETS_OPTS = --seed 0 \
	      --het_threshold $(FACETS_HET_THRESHOLD) \
	      --min_nhet $(FACETS_MIN_NHET) \
	      --snp_nbhd $(FACETS_SNP_NBHD) \
	      --ndepth_max $(FACETS_NDEPTH_MAX) \
	      $(if $(facets_pre_cval.$1),--pre_cval $(facets_pre_cval.$1),--pre_cval $(FACETS_PRE_CVAL)) \
	      $(if $(facets_cval1.$1),--cval1 $(facets_cval1.$1),--cval1 $(FACETS_CVAL1)) \
	      $(if $(facets_cval2.$1),--cval2 $(facets_cval2.$1),--cval2 $(FACETS_CVAL2)) \
	      $(if $(facets_diplogr.$1),--diplogr $(facets_diplogr.$1)) \
	      $(if $(facets_purity.$1),--purity $(facets_purity.$1))
