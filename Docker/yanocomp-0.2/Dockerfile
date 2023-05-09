FROM ubuntu:20.04

MAINTAINER Mattia Furlan <mattia.furlan@iit.it>

RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get update -qq \
    && apt-get install -y \
    build-essential \
    wget

RUN wget https://github.com/bartongroup/yanocomp/archive/refs/tags/v0.2.tar.gz \
    && tar -xvzf v0.2.tar.gz \
    && rm v0.2.tar.gz

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh \
    && chmod 777 Miniconda3-py38_4.10.3-Linux-x86_64.sh \
    && ./Miniconda3-py38_4.10.3-Linux-x86_64.sh -b -p /bin/miniconda3 \
    && /bin/miniconda3/bin/conda env create -f /yanocomp-0.2/yanocomp.yml