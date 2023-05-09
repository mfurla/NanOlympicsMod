FROM ubuntu:22.04

MAINTAINER Mattia Furlan <mattia.furlan@iit.it>

RUN apt-get update

RUN export DEBIAN_FRONTEND=noninteractive \
	&& apt-get install -y python3-pip wget

RUN wget https://github.com/bartongroup/differr_nanopore_DRS/archive/refs/tags/0.2.tar.gz \
	&& tar -xvzf 0.2.tar.gz \
	&& rm 0.2.tar.gz

RUN cd differr_nanopore_DRS-0.2 \
	&& python3 setup.py build \
	&& python3 setup.py install