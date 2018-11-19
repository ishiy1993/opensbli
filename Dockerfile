FROM python:2

ENV OPS_COMPILER=gnu
ENV OPS_INSTALL_PATH=/work/OPS/ops
ENV HDF5_INSTALL_PATH=/usr/

RUN apt-get update && \
    apt-get install -y libhdf5-dev && \
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
