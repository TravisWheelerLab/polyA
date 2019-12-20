FROM python:3.7

VOLUME [ "/code" ]
RUN pip install pipenv

WORKDIR /code
ENTRYPOINT [ "make" ]
