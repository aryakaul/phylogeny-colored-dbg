.PHONY: all cluster help conda clean cleanall cleanallall reports viewconf test format checkformat edit

SHELL := /usr/bin/env bash -eo pipefail
MAKEFLAGS += --warn-undefined-variables

.SECONDARY:
.SUFFIXES:

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!! WARNING: !! TOPDIR changes automatically to .. when run from .test/ !!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TOPDIR := $(shell if [ -d ".test" ]; then echo . ; else echo .. ; fi)

# When invoked from inside .test we need to point Snakemake back to the main Snakefile.
ifeq ($(strip $(TOPDIR)),..)
	SNAKEMAKE_PARAM_DIR := --snakefile ../workflow/Snakefile --show-failed-logs
else
	SNAKEMAKE_PARAM_DIR :=
endif

CONFIG_FILE    := config.yaml
CONDA_DIR      := $(shell awk '/^conda_dir:/ {print $$2}' $(CONFIG_FILE))
USE_CONDA      := $(shell awk '/^use_conda:/ {print $$2}' $(CONFIG_FILE))
SNAKEMAKE_JOBS ?= 64

ifeq ($(strip $(CONDA_DIR)),)
  $(error 'conda_dir' not found in $(CONFIG_FILE))
endif
ifeq ($(strip $(USE_CONDA)),)
  $(error 'use_conda' not found in $(CONFIG_FILE))
endif

CONDA_DIR_ABS  := $(abspath $(TOPDIR)/$(CONDA_DIR))
ifeq ($(filter True true,$(strip $(USE_CONDA))),True)
	SNAKEMAKE_CONDA_FLAGS := --use-conda --conda-prefix="$(CONDA_DIR_ABS)"
endif

SNAKEMAKE       ?= snakemake
SNAKEMAKE_BASE_FLAGS := -p --rerun-incomplete \
	--rerun-triggers mtime \
	--rerun-triggers params \
	--rerun-triggers input \
	--rerun-triggers code \
	$(SNAKEMAKE_PARAM_DIR) $(SNAKEMAKE_CONDA_FLAGS)
SNAKEMAKE_LOCAL_FLAGS := -j $(SNAKEMAKE_JOBS) $(SNAKEMAKE_BASE_FLAGS)
SNAKEMAKE_CLUSTER_JOBS ?= $(SNAKEMAKE_JOBS)

define REQUIRE_CONDA
  $(if $(SNAKEMAKE_CONDA_FLAGS),,$(error Conda environments are disabled via use_conda in config.yaml))
endef


######################
## General commands ##
######################

all: ## Run everything
	$(SNAKEMAKE) $(SNAKEMAKE_LOCAL_FLAGS)

cluster: ## Run everything on the cluster
	$(SNAKEMAKE) --cores 9999 -j $(SNAKEMAKE_CLUSTER_JOBS) --latency-wait 60 --restart-times 0 --cluster 'sbatch -c 1 -p short --mem=30GB -t 0-10:00:00' $(SNAKEMAKE_BASE_FLAGS)

help: ## Print help messages
	@printf "$$(grep -hE '^\S*(:.*)?##' $(MAKEFILE_LIST) \
        | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\e[36m\1\\e[0m:\2/' -e 's/^\([^#]\)/    \1/g'\
        | column -c2 -t -s : )\n"

conda: ## Create the conda environments
	$(call REQUIRE_CONDA)
	$(SNAKEMAKE) --conda-create-envs-only -d .test $(SNAKEMAKE_BASE_FLAGS)

clean: ## Clean all output archives and intermediate files
	rm -fvr output/* intermediate/* || true
	@if [ -d ".test" ]; then \
		$(MAKE) -C .test clean; \
	fi

cleanall: clean ## Clean everything but Conda, Snakemake, and input files
	rm -fvr intermediate/* || true
	@if [ -d ".test" ]; then \
		$(MAKE) -C .test cleanall; \
	fi

cleanallall: cleanall ## Clean completely everything
	rm -fvr input/* || true
	@if [ -d "$(CONDA_DIR_ABS)" ]; then rm -fvr "$(CONDA_DIR_ABS)"/*; fi
	rm -fr .snakemake/ || true
	@if [ -d ".test" ]; then \
		$(MAKE) -C .test cleanallall; \
	fi


###############
## Reporting ##
###############

viewconf: ## View configuration without comments
	@cat config.yaml \
		| perl -pe 's/ *#.*//g' \
		| grep --color='auto' -E '.*\:'
	@#| grep -Ev ^$$

reports: ## Create html report
	$(SNAKEMAKE) $(SNAKEMAKE_LOCAL_FLAGS) --report report.html
	@if [ -d ".test" ]; then \
		$(MAKE) -C .test reports; \
	fi


####################
## For developers ##
####################

test: ## Run the workflow on test data
	$(SNAKEMAKE) -d .test $(SNAKEMAKE_LOCAL_FLAGS) --show-failed-logs


format: ## Reformat all source code
	snakefmt workflow
	yapf -i --recursive workflow

checkformat: ## Check source code format
	snakefmt --check workflow
	yapf --diff --recursive workflow

edit:
	nvim -p workflow/Snakefile workflow/rules/*.smk
