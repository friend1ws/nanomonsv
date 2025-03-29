FROM ubuntu:20.04
MAINTAINER Yuichi Shiraishi <friend1ws@gmail.com> 

ENV TZ=Asia/Tokyo
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y \
    git \
    wget \
    bzip2 \
    make \
    cmake \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    python3 \
    python3-pip

RUN wget https://github.com/samtools/htslib/releases/download/1.15/htslib-1.15.tar.bz2 && \
    tar jxvf htslib-1.15.tar.bz2 && \
    cd htslib-1.15 && \
    ./configure && \
    make && \
    make install 

RUN wget http://ftp.debian.org/debian/pool/main/m/mafft/mafft_7.407-2_amd64.deb && \
    dpkg -i mafft_7.407-2_amd64.deb

RUN wget https://github.com/lbcb-sci/racon/archive/refs/tags/1.5.0.tar.gz && \
    tar zxvf racon-v1.5.0.tar.gz && \
    cd racon-v1.5.0 && mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make && make install

RUN pip3 install --upgrade setuptools

RUN pip3 install pysam==0.19.1
RUN pip3 install numpy==1.23.0
RUN pip3 install parasail==1.2.4
RUN pip3 install h5py==3.13.0

RUN wget https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/archive/v1.1.tar.gz && \
    tar zxvf v1.1.tar.gz && \
    cd Complete-Striped-Smith-Waterman-Library-1.1/src && \
    gcc -Wall -O3 -pipe -fPIC -shared -rdynamic -o libssw.so ssw.c ssw.h

ENV LD_LIBRARY_PATH /Complete-Striped-Smith-Waterman-Library-1.1/src:$LD_LIBRARY_PATH

RUN wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24.tar.bz2 && \
    tar jxvf minimap2-2.24.tar.bz2 && \
    cd minimap2-2.24 && \
    make

ENV PATH $PATH:/minimap2-2.24

RUN git clone https://github.com/friend1ws/nanomonsv.git && \
    cd nanomonsv && \
    python3 -m pip install . 


