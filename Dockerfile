FROM ubuntu:22.04
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

RUN wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    tar jxvf htslib-1.21.tar.bz2 && \
    cd htslib-1.21 && \
    ./configure && \
    make && \
    make install 

RUN wget http://ftp.debian.org/debian/pool/main/m/mafft/mafft_7.407-2_amd64.deb && \
    dpkg -i mafft_7.407-2_amd64.deb


ENV CFLAGS="-O3 -march=x86-64-v3 -mtune=generic -mno-avx512f"
ENV CXXFLAGS="-O3 -march=x86-64-v3 -mtune=generic -mno-avx512f"

RUN wget https://github.com/lbcb-sci/racon/archive/refs/tags/1.5.0.tar.gz && \
    tar zxvf 1.5.0.tar.gz && \
    cd racon-1.5.0 && mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. \
		-DCMAKE_C_FLAGS="${CFLAGS}" \ 
		-DCMAKE_CXX_FLAGS="${CXXFLAGS}" .. && \
    make && make install


RUN pip3 install --upgrade pip setuptools wheel

RUN pip3 install pysam==0.23.0
RUN pip3 install numpy==2.2.4
RUN pip3 install parasail==1.3.4
RUN pip3 install h5py==3.13.0

RUN wget https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/archive/v1.1.tar.gz && \
    tar zxvf v1.1.tar.gz && \
    cd Complete-Striped-Smith-Waterman-Library-1.1/src && \
    gcc -Wall -O3 -pipe -fPIC -shared -rdynamic -o libssw.so ssw.c ssw.h

ENV LD_LIBRARY_PATH /Complete-Striped-Smith-Waterman-Library-1.1/src:$LD_LIBRARY_PATH

RUN wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28.tar.bz2 && \
    tar jxvf minimap2-2.28.tar.bz2 && \
    cd minimap2-2.28 && \
    make

ENV PATH $PATH:/minimap2-2.28

RUN git clone https://github.com/friend1ws/nanomonsv.git && \
    cd nanomonsv && \
    pip install . && \
	pip show nanomonsv


