PR=poetry run
SRC=sp
TST=test


.PHONY: format
format:
	$(PR) black $(SRC)
	$(PR) isort $(SRC)
	$(PR) black $(TST)
	$(PR) isort $(TST)

.PHONY: check-format
check-format:
	$(PR) black --check --quiet $(SRC)
	$(PR) isort --check $(SRC)
	$(PR) black --check --quiet $(TST)
	$(PR) isort --check $(TST)

.PHONY: check-style
check-style:
	$(PR) flake8 $(SRC);
	$(PR) flake8 $(TST);

.PHONY: check-types
check-types:
	$(PR) mypy $(SRC)
	$(PR) mypy $(TST)

.PHONY: test
test:
	$(PR) pytest

.PHONY: dev-server
dev-server:
	$(PR) uvicorn sp:app --reload

.PHONY: prod-server
prod-server:
	$(PR) gunicorn sp:app \
		--workers 1 \
		--worker-class uvicorn.workers.UvicornWorker \
		--bind 0.0.0.0:8000 \
		--log-config logging.conf
