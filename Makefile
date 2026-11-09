.PHONY: install verify demo reproduce
install:
	python -m pip install -r requirements.txt
verify:
	python scripts/verify_repository.py
demo:
	python run_all.py --demo --output-dir run_outputs/demo
reproduce:
	python run_all.py --output-dir run_outputs/full
