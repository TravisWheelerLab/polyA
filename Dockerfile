FROM python:3.7

VOLUME [ "/code" ]
RUN pip install poetry

WORKDIR /code
ENTRYPOINT [ "make" ]
