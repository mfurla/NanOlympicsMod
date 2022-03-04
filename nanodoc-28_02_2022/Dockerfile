FROM ubuntu:18.04

MAINTAINER Mattia Furlan <mattia.furlan@iit.it>

RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get update -qq \
    && apt-get install -y \
    build-essential \
    wget \
    unzip \
    bzip2 \
    git \
    libidn11* \
    nano \
    less \
    bc \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.10.3-Linux-x86_64.sh \
    && chmod 777 Miniconda3-py37_4.10.3-Linux-x86_64.sh \
    && ./Miniconda3-py37_4.10.3-Linux-x86_64.sh -b -p /bin/miniconda3

RUN /bin/miniconda3/bin/conda install -c bioconda -c anaconda -c pytorch pip faiss-cpu python=3.7 numpy=1.18

RUN git clone https://github.com/uedaLabR/nanoDoc.git

RUN sed -i 's/Bio==0.0.6/Bio==1.3.3/' /nanoDoc/src/requirements.txt
RUN sed -i 's/biopython==1.72/biopython==1.79/' /nanoDoc/src/requirements.txt
RUN sed -i 's/faiss-gpu==1.5.3//' /nanoDoc/src/requirements.txt
RUN sed -i 's/numpy==1.16.2//' /nanoDoc/src/requirements.txt
RUN sed -i 's/tqdm==4.31.1//' /nanoDoc/src/requirements.txt

RUN /bin/miniconda3/bin/pip install -r /nanoDoc/src/requirements.txt

RUN mkdir /nanoDoc/weight5mer && mv /nanoDoc/weight5mer_1/* /nanoDoc/weight5mer && mv /nanoDoc/weight5mer_2/* /nanoDoc/weight5mer && rm -r /nanoDoc/weight5mer_*

RUN sed -i 's/nanoDocAnalysis.modCall(wight,param, ref, refraw,tgraw, output, chrom, chrom, start, end, strand, minreadlen)/\n    try:\n        nanoDocAnalysis.modCall(wight,param, ref, refraw,tgraw, output, chrom, chrom, start, end, strand, minreadlen)\n    except OSError:\n        pass/' /nanoDoc/src/nanoDoc.py

