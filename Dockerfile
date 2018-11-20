FROM python:2

RUN apt-get update && \
    apt-get install -y mpich libhdf5-mpich-dev && \
    export OPS_COMPILER=gnu && \
    export OPS_INSTALL_PATH=/work/OPS/ops && \
    export HDF5_INSTALL_PATH=/usr/ && \
    export MPI_INSTALL_PATH=/usr/ && \
    mkdir /work && cd /work && \
    git clone https://github.com/OP-DSL/OPS.git && \
    cd $OPS_INSTALL_PATH/c && \
    sed -i '23d' Makefile && \
    sed -i '22d' Makefile && \
    sed -i '36s/mpi//' Makefile && \
    make && \
    cd /work && \
    git clone https://github.com/opensbli/opensbli.git && \
    cd opensbli && \
    pip install -r requirements.txt && \
    make install && \
    make test

ENV OPS_COMPILER=gnu
ENV OPS_INSTALL_PATH=/work/OPS/ops
ENV HDF5_INSTALL_PATH=/usr/

