env:
	conda env remove -n interactive-cheminformatics
	conda env create -f environment.yml

dev-env:
	make env
	conda run -n interactive-cheminformatics pip install isort black flake8 mypy

format:
	isort .
	black --line-length 79 .
	flake8 --max-line-length=80 .
	mypy .
