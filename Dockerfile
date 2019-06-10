FROM ubuntu:18.04
MAINTAINER acrinklaw@lji.org

ENV PACKAGES python3-dev python3-pip hmmer ncbi-blast+ git

RUN apt-get update && \
    apt-get install -y ${PACKAGES} && \
    apt-get clean

RUN pip install bio-arc
