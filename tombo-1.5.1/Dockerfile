FROM ubuntu:18.04

MAINTAINER Mattia Furlan <mattia.furlan@iit.it>

RUN apt-get update

RUN apt-get install -y wget

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.10.3-Linux-x86_64.sh

RUN chmod 777 Miniconda3-py37_4.10.3-Linux-x86_64.sh

RUN ./Miniconda3-py37_4.10.3-Linux-x86_64.sh -b -p /bin/miniconda3

RUN /bin/miniconda3/bin/conda install --yes -c anaconda numpy=1.11.3

RUN /bin/miniconda3/bin/conda install --yes -c bioconda ont-tombo