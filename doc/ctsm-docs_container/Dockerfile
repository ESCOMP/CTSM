###################################################################
# A base linux install + packages for building CTSM documentation #
###################################################################

FROM python:3.13.2-alpine

# Update the default packages and install some other necessary ones
RUN apk update
RUN apk add --no-cache git sudo
RUN apk add --no-cache git-lfs
RUN apk add --no-cache make
ADD requirements.txt requirements.txt
RUN pip3 install --break-system-packages -r requirements.txt

# Add user etc.
RUN addgroup escomp
RUN adduser -h /home/user -G escomp -s /bin/bash -D user
RUN echo 'export USER=$(whoami)' >> /etc/profile.d/escomp.sh
RUN echo 'export PS1="[\u@escomp \W]\$ "' >> /etc/profile.d/escomp.sh
RUN echo 'user ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers.d/escomp

ENV SHELL=/bin/bash \
    LANG=C.UTF-8  \
    LC_ALL=C.UTF-8

USER user
WORKDIR /home/user
CMD ["/bin/bash", "-l"]

LABEL org.opencontainers.image.title="Container for building CTSM documentation"
LABEL org.opencontainers.image.source=https://github.com/ESCOMP/CTSM
LABEL org.opencontainers.image.version="v1.0.2e"
