Before running the pipeline:
	1) Install Snakemake: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
	2) edit line 1 of Snakefile so that each bed/bim/fam group to be run is included (comma separated)

To run the pipeline from start to finish:
	snakemake -cN		# N = number of cores available to snakemake

To force snakemake to rerun a pipeline:
	snakemake -f

To run a segment of the pipeline (and it's dependancies if needed):
	snakemake -R myrule
	
To perform a dry-run:
	snakemake -n

To draw the directed acyclic graph of jobs:
	snakemake --dag | dot -Tsvg > dag.svg
