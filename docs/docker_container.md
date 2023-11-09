# Docker container

A docker container is available, which offers ready-to-use environment for the SEIS_FILO package.  There is no need to build the software by yourself any more. 

1.  Clone Github repository.

    `git clone --depth 1 https://github.com/akuhara/SEIS_FILO.git`

2. Go inside the sofware directory.

    `cd SEIS_FILO`

3. Pull docker container.

    `docker pull akuhara/seis-filo`

4. Run docker container

    `docker run -it -v (Absolute path to directory in host machine):/wrk akuhara/seis-filo`

!!! Note    
    * With the -v option in the above procedure 4, a work direcory on a host PC is mounted on `/wrk` directory in the container.