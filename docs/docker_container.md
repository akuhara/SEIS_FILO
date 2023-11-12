# Docker Image

The SEIS_FILO package is conveniently packaged as a Docker image, readily available on Docker Hub, providing an instantly deployable environment.

1. Pull the Docker Image from Docker Hub:

    `docker pull akuhara/seis-filo`

2. Run the docker container

    `docker run -it -v (Absolute path to directory on host machine):/home/seismologist/wrk akuhara/seis-filo`

!!! Note
    * Utilize the -v option to mount a working directory from your host PC onto the /wrk directory within the container. This facilitates seamless data file integration.
    * Inside the container, your username is set as 'seismologist.'
