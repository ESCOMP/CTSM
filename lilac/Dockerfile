FROM centos:latest
LABEL description="LILAC development environment"

RUN yum install -y curl; yum upgrade -y; yum update -y;  yum clean all
RUN yum -y install wget bzip2 gcc gcc-c++ gcc-gfortran mpich make git

WORKDIR /usr/src/lilac/

RUN mkdir -p external
RUN mkdir -p ci

COPY external/esmf external/esmf
COPY external/pfunit external/pfunit
COPY ci/* ci/

# Install some remaining dependencies
ENV PATH /usr/local/miniconda/bin:$PATH
RUN ./ci/install_python.sh

# Install ESMF
# TODO: what's up with the .../lib/lib0/... maybe move this somewhere more logical?
RUN pwd
RUN ls $PWD
RUN ./ci/install_esmf.sh
ENV ESMF_CONFIG_FILE /usr/lib/libO/Linux.gfortran.64.mpiuni.default/esmf.mk

# # Install PFUNIT
RUN ./ci/install_pfunit.sh
ENV PFUNIT_INSTALL /usr/pfunit
