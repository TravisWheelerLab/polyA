FROM python:3.8

RUN pip install pipenv
COPY Pipfile .
COPY Pipfile.lock .
RUN pipenv sync --dev

VOLUME [ "/code" ]
WORKDIR /code

ENTRYPOINT [ "make" ]
