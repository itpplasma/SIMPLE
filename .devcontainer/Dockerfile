# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.217.4/containers/debian/.devcontainer/base.Dockerfile

# [Choice] Debian version (use bullseye on local arm64/Apple Silicon): bullseye, buster
# ARG VARIANT="buster"
# FROM mcr.microsoft.com/vscode/devcontainers/base:0-${VARIANT}

FROM mcr.microsoft.com/devcontainers/base:ubuntu-22.04

# ** [Optional] Uncomment this section to install additional packages. **
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
     && apt-get -y install --no-install-recommends make python3 python3-dev \
     python3-pip gfortran cmake libnetcdf-dev libnetcdff-dev libopenmpi-dev \
     libopenblas-dev wget

RUN pip install numpy scipy matplotlib f90wrap fortls

RUN cd /opt ; wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.0-beta2-linux-x86_64.tar.gz && tar -xzvf julia-1.10.0-beta2-linux-x86_64.tar.gz && ln -s `pwd`/julia-1.10.0-beta2/bin/julia /usr/local/bin/julia ; cd -
