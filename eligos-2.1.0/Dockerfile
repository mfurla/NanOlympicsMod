FROM ubuntu:18.04

MAINTAINER Mattia Furlan <mattia.furlan@iit.it>

RUN apt-get update

RUN export DEBIAN_FRONTEND=noninteractive \
	&& apt-get install -y wget

RUN wget https://gitlab.com/piroonj/eligos2/-/archive/v2.1.0/eligos2-v2.1.0.tar.gz \
	&& tar -xvzf eligos2-v2.1.0.tar.gz \
	&& rm eligos2-v2.1.0.tar.gz

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.10.3-Linux-x86_64.sh \
	&& chmod 777 Miniconda3-py37_4.10.3-Linux-x86_64.sh \
	&& ./Miniconda3-py37_4.10.3-Linux-x86_64.sh -b -p /bin/miniconda3

RUN /bin/miniconda3/bin/conda install --yes -c bioconda -c conda-forge -c anaconda python=3.6 pysam=0.13 pandas=0.23.4 pybedtools=0.8.0 bedtools=2.25 rpy2=2.8.5 r-base=3.4.1 tqdm=4.40.2 numpy=1.11.3

RUN /bin/miniconda3/bin/Rscript -e 'install.packages("samplesizeCMH", repos="https://cloud.r-project.org")'

ENV PATH=$PATH:/bin/miniconda3/bin/