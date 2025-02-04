FROM ubuntu:20.04
ENV TZ=UTC
ENV DEBIAN_FRONTEND=noninteractive

RUN apt update 
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    git \
    build-essential \
    cmake \
    wget \
    curl \
    zip \
    unzip \
    tar \
    protobuf-compiler \
    libboost-all-dev \
    pkg-config \
    capnproto 

WORKDIR /HOME

# RUN git clone https://github.com/TurakhiaLab/panman.git
# RUN panman/install/installationUbuntu.sh
