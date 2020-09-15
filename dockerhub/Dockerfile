FROM geodels/gospl-base:petsc-3.13.0
MAINTAINER Tristan Salles

ENV PETSC_DIR=/opt/petsc
ENV PETSC_ARCH=arch-linux-c-opt

RUN apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
        libgl1-mesa-dev

# install extras in a new layer
RUN python3 -m pip install --no-cache-dir \
        pytest-mpi \
        pytest-cov \
        coveralls \
        twine \
        https://github.com/j08lue/pycpt/archive/master.zip \
        pyevtk

RUN python3 -m pip install --no-cache-dir \
        itkwidgets==0.25.2 \
        pyvista==0.26.0

RUN pip3 install netCDF4
RUN apt-get update -qq && \
    apt-get install -yq --no-install-recommends \
        build-essential \
        gdal-bin \
        libgdal-dev \
        python3-gdal

# install gopspl
WORKDIR /live/lib
RUN git clone https://github.com/Geodels/gospl.git && \
    cd gospl && \
    export F90=gfortran && \
    export PETSC_DIR=/opt/petsc && \
    export PETSC_ARCH=arch-linux-c-opt && \
    python3 setup.py install

# install JIGSAW
RUN git clone https://github.com/dengwirda/jigsaw-python.git  && \
    cd jigsaw-python && \
    python3 setup.py build_external && \
    python3 setup.py install

# setup space for working in
VOLUME /live/share
WORKDIR /live

ENV LD_LIBRARY_PATH "$LD_LIBRARY_PATH:/live/lib/:/live/share"

EXPOSE 8888
ENV TINI_VERSION v0.18.0
ENTRYPOINT ["/usr/local/bin/tini", "--"]

EXPOSE 9999

# note we use xvfb which to mimic the X display
ENTRYPOINT ["/usr/local/bin/xvfbrun.sh"]

# launch notebook
ADD run-jupyter.sh /usr/local/bin/run-jupyter.sh
RUN chmod +x /usr/local/bin/run-jupyter.sh
ADD bashrc-term /root/.bashrc
#CMD /usr/local/bin/run-jupyter.sh
CMD /usr/local/bin/run-jupyterlab.sh
