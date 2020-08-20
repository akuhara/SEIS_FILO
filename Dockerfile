From ubuntu:18.04
ARG USER=seismologist
RUN apt-get update && apt-get install -y \
  gfortran \
  liblapack-dev \
  libopenmpi-dev \
  libfftw3-dev \
  openmpi-bin \
  python3-pip \
  ssh
COPY requirements.txt /tmp/
RUN pip3 install --upgrade pip && pip3 install --requirement /tmp/requirements.txt
COPY . /usr/local/SEIS_FILO/
WORKDIR /usr/local/SEIS_FILO/src
RUN make FFTW="-I/usr/include -lfftw3"
RUN useradd -m ${USER}
USER ${USER}
WORKDIR /home/${USER}
ENV PATH $PATH:/usr/local/SEIS_FILO/bin
RUN cp -r /usr/local

CMD ["/bin/bash"] 

