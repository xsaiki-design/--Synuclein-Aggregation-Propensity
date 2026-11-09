.PHONY: help setup data data-synth inspect analysis test check all clean
PY := PYTHONPATH=src python

help:
	@echo "setup       install dependencies (pip install -e .[dev])"
	@echo "data        build the REAL curated CSV from data/raw/SI.pdf"
	@echo "data-synth  write the synthetic stand-in so the pipeline runs now"
	@echo "inspect     dump the SI.pdf tables for review (needs data/raw/SI.pdf)"
	@echo "analysis    run decomposition + prediction, write figures + results.md"
	@echo "test        run the test suite"
	@echo "check       test + analysis (one-command verify)"
	@echo "all         data-synth + analysis"

setup:
	pip install -e ".[dev]"

data:
	$(PY) -m asyn.ingest

data-synth:
	$(PY) -m asyn.synth

inspect:
	$(PY) -m asyn.ingest --inspect

analysis:
	$(PY) -m asyn.pipeline

test:
	$(PY) -m pytest -q

check: test analysis

all: data-synth analysis

clean:
	rm -rf reports/figures/*.png reports/results.md data/raw/si_tables_dump.csv
