FROM ubuntu:22.04

MAINTAINER Mattia Furlan <mattia.furlan@iit.it>

RUN export DEBIAN_FRONTEND=noninteractive \
	&& apt-get update && \
	apt-get -y install \
		autoconf \
		python3 \
		python3-pip \
		python3-numpy \
		python3-h5py \
		python3-six \
		python3-progressbar

RUN pip install ont-fast5-api
