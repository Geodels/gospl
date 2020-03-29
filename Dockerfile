FROM geodels/gospl-base:latest

ENV PETSC_DIR=/opt/petsc
ENV PETSC_ARCH=arch-linux-c-opt
ENV LD_LIBRARY_PATH "$LD_LIBRARY_PATH:/live/lib/"
ADD dockerhub/bashrc-term /root/.bashrc

# install gopspl
WORKDIR /live/lib
RUN git clone https://github.com/Geodels/gospl.git && \
    cd gospl && \
    export F90=gfortran && \
    export PETSC_DIR=/opt/petsc && \
    export PETSC_ARCH=arch-linux-c-opt && \
    python3 setup.py install

WORKDIR /live/lib
RUN pip3 install -e .

CMD ["jupyter", "notebook", "--ip='0.0.0.0'", "--no-browser"]
