FROM python:3.8

VOLUME [ "/code" ]
RUN pip install poetry

WORKDIR /code
ENTRYPOINT [ "make" ]
