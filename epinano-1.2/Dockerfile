FROM ubuntu:18.04

MAINTAINER Mattia Furlan <mattia.furlan@iit.it>

RUN apt-get update
RUN apt update

RUN export DEBIAN_FRONTEND=noninteractive \
	&& apt-get install -y wget python3-pip python3-venv gcc git default-jdk samtools curl libxml2-dev libssl-dev libcurl4-openssl-dev

RUN pip3 install --upgrade pip
RUN pip3 install atomicwrites==1.4.0 attrs==21.2.0 biopython==1.76 cloudpickle==1.3.0
RUN pip3 install dask==2.5.2 fsspec==2021.6.1 future==0.17.1 h5py==2.10.0 importlib-metadata==4.6.1
RUN pip3 install locket==0.2.1  more-itertools==8.8.0 numpy==1.17.2 pandas==0.24.2
RUN pip3 install partd==1.2.0 pluggy==0.13.1 py==1.10.0 pysam==0.15.4 pytest==4.4.1 python-dateutil==2.8.1
RUN pip3 install pytz==2021.1 scikit-learn==0.20.2 scipy==1.5.4 six==1.16.0 toolz==0.11.1 typing-extensions==3.10.0.0 zipp==3.5.0

RUN wget https://github.com/novoalab/EpiNano/archive/refs/tags/Epinano1.2.1.tar.gz \
	&& tar -xvzf Epinano1.2.1.tar.gz \
	&& rm Epinano1.2.1.tar.gz

RUN git clone "https://github.com/lindenb/jvarkit.git" \
	&& cd jvarkit \
	&& ./gradlew sam2tsv

RUN sed -i 's/java -jar  {sam2tsv} -r {reference_file}/java -jar  {sam2tsv} -R {reference_file}/' /EpiNano-Epinano1.2.1/Epinano_Variants.py

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.10.3-Linux-x86_64.sh \
	&& chmod 777 Miniconda3-py37_4.10.3-Linux-x86_64.sh \
	&& ./Miniconda3-py37_4.10.3-Linux-x86_64.sh -b -p /bin/miniconda3 \
	&& /bin/miniconda3/bin/conda install --yes -c conda-forge -c r r-base=3.6 r-rcurl r-outliers r-reshape2 \
	r-ggplot2 r-car r-ggrepel r-tidyverse r-optparse