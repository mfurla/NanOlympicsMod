FROM ubuntu:22.04

MAINTAINER Mattia Furlan <mattia.furlan@iit.it>

RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get update -qq \
    && apt-get install -y \
    build-essential \
    wget \
    unzip \
    bzip2 \
    git \
    libidn11* \
    python3-pip \
    samtools \
    bedtools \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install seaborn scipy pandas numpy biopython matplotlib statsmodels
RUN git clone https://github.com/DepledgeLab/DRUMMER.git