run:
	python3 scripts/run_experiments.py

sweep:
	python3 numerics/simulations/parameter_sweep.py

periodic_twisted_l1:
	python3 numerics/simulations/periodic_twisted_l1.py

periodic_twisted_l1_flux_scan:
	python3 numerics/simulations/periodic_twisted_l1_flux_scan.py

periodic_twisted_l1_hodge_projection:
	python3 numerics/simulations/periodic_twisted_l1_hodge_projection.py
