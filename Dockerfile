
FROM python:3.7-slim-bullseye AS easel

ARG easel_version=0.48

RUN apt-get update && apt-get -y install \
    autoconf \
    build-essential
ADD https://github.com/EddyRivasLab/easel/archive/easel-${easel_version}.tar.gz easel.tar.gz
RUN tar xf easel.tar.gz
WORKDIR /easel-easel-${easel_version}
RUN autoconf && \
    ./configure && \
    make && \
    gcc -g -Wall -I. -L. -o esl_scorematrix \
        -DeslSCOREMATRIX_EXAMPLE \
        esl_scorematrix.c -leasel -lm && \
    cp esl_scorematrix /

# -----------------------------------------------------------------------------

FROM python:3.7-slim-bullseye as ultra

RUN apt-get update && apt-get -y install \
    git \
    build-essential \
    cmake
RUN git clone https://github.com/TravisWheelerLab/ULTRA.git
WORKDIR /ULTRA
RUN cmake . && \
    make && \
    cp ultra /

# -----------------------------------------------------------------------------

FROM python:3.7-slim-bullseye

# Install and test esl_scorematrix
COPY --from=easel /esl_scorematrix /usr/local/bin/
RUN esl_scorematrix -h

# Install and test ultra
COPY --from=ultra /ultra /usr/local/bin/
RUN ultra -h

COPY requirements.txt .
RUN pip install -r requirements.txt

VOLUME [ "/code" ]
WORKDIR /code

ENTRYPOINT [ "make" ]
