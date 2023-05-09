FROM ubuntu:18.04

MAINTAINER Mattia Furlan <mattia.furlan@iit.it>

RUN export DEBIAN_FRONTEND=noninteractive \
        && apt-get update && \
        apt-get -y install \
        wget \
        git
        
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.10.3-Linux-x86_64.sh

RUN chmod 777 Miniconda3-py37_4.10.3-Linux-x86_64.sh

RUN ./Miniconda3-py37_4.10.3-Linux-x86_64.sh -b -p /bin/miniconda3

RUN git clone https://github.com/YeoLab/MINES.git

RUN sed -i 's/1.1.1=h7b6447c_0/1.1.1/' /MINES/environment.yml

RUN /bin/miniconda3/bin/conda env create -f /MINES/environment.yml

ENV PATH=$PATH:/bin/miniconda3/envs/mines/bin/