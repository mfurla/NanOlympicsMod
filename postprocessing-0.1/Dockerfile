FROM ubuntu:22.04

MAINTAINER Mattia Furlan <mattia.furlan@iit.it>

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq && \        
	apt-get install -y \
	r-base \
	libcurl4-openssl-dev \
	libxml2-dev \
	libssl-dev

RUN R -e "install.packages('xml2')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('IRanges')"  
RUN R -e "BiocManager::install('GenomicRanges')" 
RUN R -e "BiocManager::install('ensembldb')"