FROM geodels/gospl-base:latest

ENV PETSC_DIR=/opt/petsc
ENV PETSC_ARCH=arch-linux-c-opt
ENV LD_LIBRARY_PATH "$LD_LIBRARY_PATH:/live/lib/"
ADD dockerhub/bashrc-term /root/.bashrc

# VTK openGL driver
RUN apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
        libgl1-mesa-dev

# install extras in a new layer
RUN python3 -m pip install --no-cache-dir \
        pytest-mpi \
        pytest-cov \
        coveralls

# install gopspl
WORKDIR /live/lib

RUN git clone https://github.com/Geodels/gospl.git && \
    cd gospl && \
    export F90=gfortran && \
    export PETSC_DIR=/opt/petsc && \
    export PETSC_ARCH=arch-linux-c-opt && \
    python3 setup.py install && \
    python3 -m pip install -e .

CMD ["jupyter", "notebook", "--ip='0.0.0.0'", "--no-browser"]
