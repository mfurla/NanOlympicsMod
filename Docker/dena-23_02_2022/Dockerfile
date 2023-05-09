FROM ubuntu:20.04

RUN export DEBIAN_FRONTEND=noninteractive \
        && apt-get update && \
        apt-get -y install \
                autoconf \
                python3 \
                python3-pip \
                python3-tk \
                python3-numpy \
                python3-scipy \
                python3-pandas \
                git \
                wget \
                lsb-release \
                apt-transport-https \
                vim \
                zlib1g-dev \
                libbz2-dev \
                liblzma-dev \
                libcurl4-openssl-dev \
                libidn11 \
                curl \
                unzip

RUN pip3 install \
    sklearn \
    pysam \
    torch \
    torchvision \
    ont-tombo

RUN git clone https://github.com/weir12/DENA.git

RUN wget https://qpoviq.DM.files.1drv.com/y4mUf8TBx4Ji7PjPRvldPKYBJyfbZSXxo7beR6J-m98d3oPYK3YVr_Ah6zFwtdeq_06pZ0uvCMH-jsRSq5yeg9FA80H-SE_uPamEiVLCDLMe4fqKkVz8z1eMLckE7yvrnLHdf37qAsMInIQVcussqVDw_nxoMCjh1giHfc7lW6HL76cG2c4bW6OViG_GWbbTJyMF3axQbkMqStR6M5Wlrh3NQ -O denaModels.zip \
&& mkdir DENA/denaModels/ \
&& mv denaModels.zip DENA/denaModels \
&& cd DENA/denaModels \
&& unzip denaModels.zip