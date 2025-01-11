# hash:sha256:ab3d00bdb30f7374b33ac821722b97b13bbcfbbb903f42bf4398a8338b6caf25
FROM registry.codeocean.com/codeocean/miniconda3:4.12.0-python3.9-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN pip3 install -U --no-cache-dir \
    biopython==1.83 \
    more-itertools==10.2.0 \
    openpyxl==3.1.2 \
    pandas==2.2.0 \
    scipy==1.12.0 \
    statsmodels==0.14.1 \
    swalign==0.3.7
