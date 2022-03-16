FROM ubuntu:18.04 as first_base

LABEL mantainer="Livia Moura"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y \
  curl \
  g++ \
  libbz2-dev \
  liblzma-dev \
  make \
  python \
  tar \
  tcllib \
  unzip \
  wget \
  zlib1g-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libncurses5-dev \
  git

# samtools 
ADD https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 samtools-1.14.tar.bz2
RUN tar -xf samtools-1.14.tar.bz2 &&\
    cd samtools-1.14 && \
    make

# bcftools
ADD https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2 bcftools-1.14.tar.bz2
RUN tar -xf bcftools-1.14.tar.bz2 &&\
    cd bcftools-1.14 &&\
    make

# bbmap

ADD https://ufpr.dl.sourceforge.net/project/bbmap/BBMap_38.96.tar.gz BBMap_38.96.tar.gz
RUN tar -xf BBMap_38.96.tar.gz

# FixAME
RUN git clone https://github.com/LiviaMoura/FixAME.git

### Final repository
FROM openjdk:19-slim as final_repository
COPY --from=python:3.6.15-slim / /

LABEL mantainer="Livia Moura"

RUN apt-get update &&\
    apt-get -y install vim libcurl4-openssl-dev libncurses5 &&\
    apt-get clean &&\
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip install --upgrade pip && pip install regex xopen==0.9.0 pandas==1.1.5 pysam==0.15.4 biopython==1.77

COPY --from=first_base /samtools-1.14/samtools  /usr/bin/
COPY --from=first_base /bcftools-1.14/bcftools  /usr/bin/
COPY --from=first_base /bbmap/ /usr/bin/
COPY --from=first_base /FixAME/ /usr/local/bin/

RUN chmod +x /usr/local/bin/FixAME.py && echo "alias FixAME=FixAME.py" >> /root/.bashrc

CMD FixAME