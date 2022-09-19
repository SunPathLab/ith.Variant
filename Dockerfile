FROM ubuntu:bionic

## Docker image for ith.variant

ENV DEBIAN_FRONTEND noninteractive

LABEL maintainer="Aziz Khan (azez.khan@gmail.com)"

ARG ZLIB_VERSION=1.2.11
ARG BAMTOOLS_VERSION=2.5.1
ARG BOOST_VERSION=1.54.0

# install required packages
RUN apt-get update \
    && apt-get install -y locales git pkg-config libjsoncpp-dev \
        build-essential \
        libz-dev \
        wget \
        make \
        cmake \
        gcc \
        g++ \
        zlib1g-dev

# Set the locale
RUN locale-gen en_US.UTF-8  
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8

WORKDIR /opt/

# install zlib
RUN cd /opt/ && \
	wget http://www.zlib.net/fossils/zlib-${ZLIB_VERSION}.tar.gz && \
	tar -xzf zlib-${ZLIB_VERSION}.tar.gz && \
	mv zlib-${ZLIB_VERSION} zlib && cd zlib && \
	./configure --prefix=/opt/zlib && \
	make install && \
	cd /opt/ && rm -rf zlib-${ZLIB_VERSION}*

ENV PKG_CONFIG_PATH /opt/zlib/lib/pkgconfig
ENV LD_LIBRARY_PATH /opt/zlib/lib:${LD_LIBRARY_PATH}

# install bamtools
RUN cd /opt/ && \
	wget --quiet https://github.com/pezmaster31/bamtools/archive/refs/tags/v${BAMTOOLS_VERSION}.tar.gz \
	-O bamtools-${BAMTOOLS_VERSION}.tar.gz && \
	tar -xzf bamtools-${BAMTOOLS_VERSION}.tar.gz && \
	rm bamtools-${BAMTOOLS_VERSION}.tar.gz \
	&& mv bamtools-* bamtools && cd bamtools \
	&& mkdir build \
	&& cd build \
	&& cmake -DCMAKE_INSTALL_PREFIX=/opt/bamtools/ .. \
	&& make \
	&& make install

# ith.Variant look for api in this dir
RUN cp /opt/bamtools/include/bamtools/* /opt/bamtools/include/ -r

## install boost
RUN cd /opt/ && \
	wget --no-check-certificate https://cfhcable.dl.sourceforge.net/project/boost/boost/1.54.0/boost_1_54_0.tar.gz \
	&& tar -xzf boost_1_54_0.tar.gz \
	&& cd /opt/boost_1_* \
	&& ./bootstrap.sh \
	&& ./b2 install --with-regex --prefix=/opt/boost \
	&& rm /opt/boost_* -rf

ENV BOOST_ROOT=/opt/boost/
ENV BAMTOOLS_ROOT=/opt/bamtools/
ENV ZLIB_ROOT=/opt/zlib/

# install ith.Variant
RUN cd /opt/ && \
	git clone https://github.com/SunPathLab/ith.Variant.git && cd ith.Variant \
	&& make BAMTOOLS_ROOT=/opt/bamtools/ ZLIB_ROOT=/opt/zlib/ BOOST_ROOT=/opt/boost/

## Copy local files
COPY requirements.txt /opt/

CMD [ "/bin/bash" ]

# install miniconda
ARG CONDA_VERSION=py37_4.10.3
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh -O ~/miniconda.sh && \
	bash ~/miniconda.sh -b -p /opt/conda && \
	rm ~/miniconda.sh

ENV PATH /opt/conda/bin:$PATH

# update conda env
RUN conda update -n base -c defaults conda -y

# install other ith.Variant requirenments
RUN conda install -c conda-forge mamba --yes \
	#&& activate ith.variant \
	&& mamba install --file /opt/requirements.txt -c bioconda -c conda-forge -c defaults -y

# symbolic link perl to conda
RUN ln -sf /opt/conda/bin/perl /usr/bin/perl 

# install TitanCNA modified version
COPY ./pkgs/TitanCNA_1.26.0.tar.gz /opt/
RUN R -e 'install.packages("/opt/TitanCNA_1.26.0.tar.gz")'
RUN rm /opt/TitanCNA_1.26.0.tar.gz

ENV PATH=/opt:/opt/ith.Variant/bin:${PATH}
