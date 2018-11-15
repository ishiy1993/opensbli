FROM python:2

RUN mkdir /work && cd /work && \
    git clone https://github.com/OP-DSL/OPS.git && \
    export OPS_COMPILER=gnu && \
    export OPS_INSTALL_PATH=/work/OPS/ops && \
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

 

