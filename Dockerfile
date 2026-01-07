FROM ubuntu:22.04
LABEL maintainer="Yuichi Shiraishi <friend1ws@gmail.com>" 

ENV TZ=Asia/Tokyo
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y \
    git \
    wget \
    bzip2 \
    make \
    cmake \
    gcc \
    g++ \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    python3 \
    python3-pip \
	mafft

RUN wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    tar jxvf htslib-1.21.tar.bz2 && \
    cd htslib-1.21 && \
    ./configure && \
    make && \
    make install 


ENV CFLAGS="-O3"
ENV CXXFLAGS="-O3"

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
RUN pip3 install boto3==1.37.28


RUN wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28.tar.bz2 && \
    tar jxvf minimap2-2.28.tar.bz2 && \
    cd minimap2-2.28 && \
    make

ENV PATH=$PATH:/minimap2-2.28

RUN git clone https://github.com/friend1ws/nanomonsv.git && \
    cd nanomonsv && \
    pip install . && \
	pip show nanomonsv


