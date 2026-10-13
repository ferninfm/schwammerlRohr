FROM ubuntu:22.04
RUN apt-get update && apt-get install -y --no-install-recommends wget ca-certificates && \
    wget -q https://github.com/bbuchfink/diamond/releases/download/v2.1.12/diamond-linux64.tar.gz -O /tmp/diamond.tar.gz && \
    tar -xzf /tmp/diamond.tar.gz -C /usr/local/bin && chmod +x /usr/local/bin/diamond && \
    rm /tmp/diamond.tar.gz && apt-get clean && rm -rf /var/lib/apt/lists/*
