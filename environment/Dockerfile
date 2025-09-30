# hash:sha256:5bb106099fb3e06a246854868a014f221f25b89f63dbd6176e67f45b1116b592
FROM registry.codeocean.com/codeocean/miniconda3:4.12.0-python3.9-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        fort77=1.15-11 \
        gcc=4:9.3.0-1ubuntu2 \
        gfortran=4:9.3.0-1ubuntu2 \
        libgfortran5=10.5.0-1ubuntu1~20.04 \
        libgmp3-dev=2:6.2.0+dfsg-4ubuntu0.1 \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y --channel=conda-forge --channel=https://conda.rosettacommons.org \
        numpy==1.26.4 \
        pandas==2.2.3 \
        pyrosetta==2025.06+release.029c6a159b \
    && conda clean -ya

COPY postInstall /
RUN /postInstall
