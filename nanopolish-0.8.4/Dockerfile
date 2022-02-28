FROM ubuntu:22.04

MAINTAINER Mattia Furlan <mattia.furlan@iit.it>

RUN export DEBIAN_FRONTEND=noninteractive \
	&& apt-get update && \
	apt-get -y install \
		git \
		build-essential \
		wget \
		libz-dev

RUN git clone --recursive https://github.com/jts/nanopolish.git \
	&& cd nanopolish \
	&& make

RUN mv /nanopolish /bin
