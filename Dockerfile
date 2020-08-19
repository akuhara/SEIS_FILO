From ubuntu:18.04
RUN apt-get update && apt-get install -y \
  gfortran \
  liblapack-dev \
  libopenmpi-dev \
  libfftw3-dev \
  openmpi-bin \
  python3-pip
COPY requirements.txt /tmp/
RUN pip install --requirement /tmp/requirements.txt
COPY src/ /
CMD ["/bin/bash"] 
  
