FROM python:3.8 AS easel

ARG easel_version=0.48

RUN apt-get update
RUN apt-get -y install \
    autoconf \
    build-essential \
    wget

RUN wget https://github.com/EddyRivasLab/easel/archive/easel-${easel_version}.tar.gz
RUN tar xf easel-${easel_version}.tar.gz
WORKDIR easel-easel-${easel_version}

RUN autoconf
RUN ./configure
RUN make
# Build the bit of Easel we need
RUN gcc -g -Wall -I. -L. -o esl_scorematrix \
    -DeslSCOREMATRIX_EXAMPLE \
    esl_scorematrix.c -leasel -lm

RUN cp esl_scorematrix /

FROM python:3.8 as ultra

RUN apt-get update
RUN apt-get -y install \
    autoconf \
    build-essential \
    wget

RUN git clone https://github.com/TravisWheelerLab/ULTRA.git
WORKDIR /ULTRA
RUN git checkout prerelease
RUN make

RUN cp ultra /

FROM python:3.8

# Install and test esl_scorematrix
COPY --from=easel /esl_scorematrix /usr/local/bin/
RUN esl_scorematrix -h

# Install and test ultra
COPY --from=ultra /ultra /usr/local/bin/
RUN ultra -h

RUN pip install pipenv
COPY Pipfile .
COPY Pipfile.lock .
RUN pipenv sync --dev

VOLUME [ "/code" ]
WORKDIR /code

ENTRYPOINT [ "make" ]
