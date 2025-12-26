FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        git \
        zlib1g-dev \
        liblzma-dev \
        libbz2-dev \
        libcurl4-openssl-dev \
        libssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/dmtools
COPY . .

RUN make -j"$(nproc)"

ENV PATH="/opt/dmtools:${PATH}"

CMD ["./dmtools", "--help"]
