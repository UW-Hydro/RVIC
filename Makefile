# Used for GDAL installation
export PIP_INDEX_URL=https://pypi.pacificclimate.org/simple

# Setup venv var
ifeq ($(TMPDIR),)
VENV_PATH := /tmp/rvic-daccs-venv
else
VENV_PATH := $(TMPDIR)/rvic-daccs-venv
endif

# Makefile Vars
SHELL:=/bin/bash
PYTHON=${VENV_PATH}/bin/python3
PIP=${VENV_PATH}/bin/pip

.PHONY: all
all: install pre-commit-hook test

.PHONY: clean-venv
clean-venv:
	rm -rf $(VENV_PATH)

.PHONY: install
install: venv
	${PIP} install -U pip
	${PIP} install -r requirements.txt -r test_requirements.txt
	${PIP} install -e .

.PHONY: install-ci
install-ci:
	pip3 install -U pip
	pip3 install -r requirements.txt -r test_requirements.txt
	pip3 install -e .

.PHONY: pre-commit-hook
pre-commit-hook: venv
	source ${VENV_PATH}/bin/activate
	pip3 install pre-commit
	pre-commit install

.PHONY: test
test: venv
	${PYTHON} -m pytest -vv

.PHONY: venv
venv:
	test -d $(VENV_PATH) || python3 -m venv $(VENV_PATH)
