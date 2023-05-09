FROM ubuntu:22.04

MAINTAINER Mattia Furlan <mattia.furlan@iit.it>

RUN apt-get update

RUN export DEBIAN_FRONTEND=noninteractive \
	&& apt-get install -y wget python3-pip git libncurses5

RUN pip3 install gdown

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
	&& chmod 777 Miniconda3-latest-Linux-x86_64.sh \
	&& ./Miniconda3-latest-Linux-x86_64.sh -b -p /bin/miniconda3

RUN /bin/miniconda3/bin/conda install --yes -c hcc jvarkit-sam2tsv

RUN gdown https://drive.google.com/uc?id=18HtAc358ks0pcInwS1RUNhWT-vh6k608

RUN tar -xvzf nanom6A_2021_10_22.tar.gz \
	&& rm nanom6A_2021_10_22.tar.gz

ENV PATH=$PATH:/bin/miniconda3/bin/

RUN git clone https://github.com/broadinstitute/picard.git \
	&& cd picard/ \
	&& ./gradlew shadowJar