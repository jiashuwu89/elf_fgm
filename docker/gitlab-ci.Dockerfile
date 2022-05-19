FROM python:3.9-slim

RUN apt-get update \
    && apt-get install --no-install-recommends -y make curl

ENV PATH="/root/.local/bin:$PATH"

RUN curl -sSL https://install.python-poetry.org | python3 - \
    && poetry config virtualenvs.in-project true
