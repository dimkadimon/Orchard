PY ?= python3

.PHONY: all bounds validate figures paper clean

all: bounds validate figures

bounds:
	$(PY) second_order_k.py
	$(PY) second_order_full.py

validate:
	$(PY) verify_second_order.py
	$(PY) validate_theorem.py
	$(PY) validate_k.py

figures:
	$(PY) make_figures.py
	$(PY) fig_t4b.py

paper:
	cd paper && tectonic -X compile orchard.tex

clean:
	rm -rf __pycache__ paper/*.aux paper/*.log paper/*.out paper/*.synctex.gz
