FROM centos:latest
LABEL description="LILAC development environment"


RUN yum install -y curl; yum upgrade -y; yum update -y;  yum clean all
RUN yum -y install wget bzip2 gcc gcc-c++ gcc-gfortran mpich make

ADD ./ /usr/src/lilac
RUN cd /usr/src/lilac && mkdir -p build

WORKDIR /usr/src/lilac
ENV PATH /usr/local/miniconda/bin:$PATH
RUN /usr/src/lilac/ci/install_python.sh
RUN /usr/src/lilac/ci/install_esmf.sh
RUN /usr/src/lilac/ci/install_pfunit.sh

# RUN mkdir -p /usr/src/lilac/build && cd /usr/src/lilac/build && cmake ..

CMD /bin/bash -c "ctest"
