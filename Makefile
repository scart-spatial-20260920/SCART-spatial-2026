# SCART spatial transcriptomics - analysis pipeline
#
#   export SCART_DATA=/path/to/geo/download
#   make all
#
# Every target reads from $(SCART_DATA) and writes to $(RESULTS).  Nothing writes
# back into the data directory.

SCART_DATA ?= ./data
RESULTS    ?= ./results
PY         ?= python

export SCART_DATA
export RESULTS

SCRIPTS := scripts

.PHONY: all check load celltype neighborhood deg pathway morans figures fig4 audit geo clean help

help:
	@echo "make check       - verify environment and that SCART_DATA exists"
	@echo "make all         - full pipeline"
	@echo "make figures     - all journal figures from cached intermediates"
	@echo "make fig4        - regenerate Figure 4 (and panel C alone)"
	@echo "make audit       - zone-count and normalisation audit tables"
	@echo "make geo         - build the GEO submission files"
	@echo "make clean       - remove $(RESULTS) (never touches SCART_DATA)"

check:
	@$(PY) -c "import scanpy,anndata,squidpy,decoupler,matplotlib; print('env ok')"
	@test -d "$(SCART_DATA)" || (echo "ERROR: SCART_DATA=$(SCART_DATA) does not exist"; exit 1)
	@echo "data ok: $(SCART_DATA)"

all: check load celltype neighborhood deg pathway morans figures audit

load:
	$(PY) $(SCRIPTS)/phase2_03_scanpy_load.py

celltype:
	$(PY) $(SCRIPTS)/phase2_04_spatial_celltype.py

neighborhood:
	$(PY) $(SCRIPTS)/phase2_05_squidpy_neighborhood.py

deg:
	$(PY) $(SCRIPTS)/phase2_06_deg_rim.py

pathway:
	$(PY) $(SCRIPTS)/phase2_07_pathway_decoupler.py

morans:
	$(PY) $(SCRIPTS)/phase2_08_morans_i.py

figures:
	$(PY) $(SCRIPTS)/phase2_10_journal_figures.py

fig4:
	$(PY) $(SCRIPTS)/regenerate_fig4C.py

audit:
	$(PY) $(SCRIPTS)/zone_truth.py
	$(PY) $(SCRIPTS)/normalisation_audit.py

geo:
	$(PY) $(SCRIPTS)/export_geo_matrices.py \
		--h5ad "$(SCART_DATA)/SCART_spatial_annotated.h5ad" \
		--out  "$(RESULTS)/SCART_GEO_submission"

clean:
	rm -rf "$(RESULTS)"
