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

RUN apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
        nodejs \
        npm \
        libgl1-mesa-dev \
        xvfb

RUN python3 -m pip install --no-cache-dir \
        voila \
        panel \
        jupyterlab

RUN jupyter labextension install @jupyter-voila/jupyterlab-preview
RUN jupyter labextension install @pyviz/jupyterlab_pyviz

# install JIGSAW
RUN git clone https://github.com/dengwirda/jigsaw-python.git  && \
cd jigsaw-python && \
python3 setup.py build_external && \
python3 setup.py install

RUN python3 -m pip install --no-cache-dir \
        itkwidgets \
        pyvista

RUN pip3 install netCDF4
RUN apt-get update -qq && \
    apt-get install -yq --no-install-recommends \
        build-essential \
        gdal-bin \
        libgdal-dev \
        python3-gdal


RUN jupyter labextension install @jupyter-widgets/jupyterlab-manager \
    jupyter-matplotlib jupyterlab-datawidgets itkwidgets

# install gopspl
WORKDIR /live/lib
RUN git clone https://github.com/Geodels/gospl.git   && \
    cd gospl && \
    export F90=gfortran && \
    export PETSC_DIR=/opt/petsc && \
    export PETSC_ARCH=arch-linux-c-opt && \
    python3 setup.py install

# setup space for working in
VOLUME /live/share
WORKDIR /live
RUN rm -rf jigsaw-python

ENV LD_LIBRARY_PATH "$LD_LIBRARY_PATH:/live/lib/:/live/share"

EXPOSE 8888
COPY start_xvfb.sh /usr/local/sbin/start_xvfb.sh
RUN chmod a+x /usr/local/sbin/start_xvfb.sh

ENTRYPOINT ["tini", "-g", "--", "start_xvfb.sh"]

EXPOSE 9999

RUN pip3 install lavavu
# note we use xvfb which to mimic the X display for lavavu
ENTRYPOINT ["/usr/local/bin/xvfbrun.sh"]


# launch jupyter
#ADD run-jupyter.sh /usr/local/bin/run-jupyter.sh
#RUN chmod +x /usr/local/bin/run-jupyter.sh
#ADD bashrc-term /root/.bashrc
#CMD /usr/local/bin/run-jupyter.sh

# launch jupyterlab
ADD run-jupyterlab.sh /usr/local/bin/run-jupyterlab.sh
RUN chmod +x /usr/local/bin/run-jupyterlab.sh
ADD bashrc-term /root/.bashrc
CMD /usr/local/bin/run-jupyterlab.sh
