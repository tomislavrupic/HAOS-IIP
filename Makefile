run:
	python3 scripts/run_experiments.py

pdf:
	python3 scripts/release_phase.py $(if $(PHASE_DIR),--phase-dir "$(PHASE_DIR)",) $(if $(PAPER),--paper "$(PAPER)",) $(if $(TITLE),--title "$(TITLE)",) $(if $(COMMIT_MSG),--commit-message "$(COMMIT_MSG)",) $(foreach path,$(EXTRA),--extra "$(path)") $(if $(DRY_RUN),--dry-run,)

pin-pdf:
	python3 scripts/release_phase.py $(if $(PHASE_DIR),--phase-dir "$(PHASE_DIR)",) $(if $(PAPER),--paper "$(PAPER)",) $(if $(TITLE),--title "$(TITLE)",) $(if $(COMMIT_MSG),--commit-message "$(COMMIT_MSG)",) $(foreach path,$(EXTRA),--extra "$(path)") --write-target

release: pdf

sweep:
	python3 numerics/simulations/parameter_sweep.py

periodic_twisted_l1:
	python3 numerics/simulations/periodic_twisted_l1.py

periodic_twisted_l1_flux_scan:
	python3 numerics/simulations/periodic_twisted_l1_flux_scan.py

periodic_twisted_l1_hodge_projection:
	python3 numerics/simulations/periodic_twisted_l1_hodge_projection.py
