format:
	isort .
	black --line-length 79 .
	flake8 --max-line-length=80 .
	mypy .
