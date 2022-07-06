FROM python:3

# Packages update
RUN apt-get update \
 && apt-get install -yq --no-install-recommends \
    openbabel \
    ca-certificates \
    build-essential \
    cmake \
    wget \
    libboost-dev \
    libboost-system-dev \
    libboost-thread-dev \
    libboost-serialization-dev \
    libboost-python-dev \
    libboost-regex-dev \
    libboost-iostreams-dev \
    libcairo2-dev \
    libeigen3-dev \
    python3-dev \
    python3-numpy \
    swig \
    python3-rdkit \
    librdkit1 \
    rdkit-data \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Build openbabel in python
ARG OPENBABEL_VERSION=openbabel-2-4-0
RUN wget --quiet https://github.com/openbabel/openbabel/archive/${OPENBABEL_VERSION}.tar.gz \
&& tar -xzf ${OPENBABEL_VERSION}.tar.gz \
&& mv openbabel-${OPENBABEL_VERSION} openbabel \
&& rm ${OPENBABEL_VERSION}.tar.gz

RUN mkdir /openbabel/build
WORKDIR /openbabel/build

RUN cmake \
  -D RUN_SWIG=ON \
  -D PYTHON_EXECUTABLE=/usr/bin/python \
  -D PYTHON_INCLUDE_DIR=/usr/include/python \
  -D CMAKE_INSTALL_PREFIX=/usr \
  -D CMAKE_BUILD_TYPE=Release \
  ..

RUN make -j $(nproc) \
 && make install

# Install python requirements
WORKDIR /basf
COPY requirements.txt /basf/
RUN pip install -r requirements.txt
COPY . /basf/