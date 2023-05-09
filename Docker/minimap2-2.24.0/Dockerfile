FROM ubuntu:22.04

MAINTAINER Mattia Furlan <mattia.furlan@iit.it>

RUN apt-get update

RUN apt-get install -y curl bzip2 samtools

RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 | tar jxf -

RUN mkdir /bin/minimap2 \
	&& mv minimap2-2.24_x64-linux/* /bin/minimap2 \
	&& rm -r minimap2-2.24_x64-linux