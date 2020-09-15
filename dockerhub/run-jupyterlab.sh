cd /live

export PATH=$PATH:/opt/mpich/bin/:/opt/petsc/bin
export SHELL=/bin/bash

jupyter lab --port=8888 --ip='*' --no-browser --allow-root --NotebookApp.token=''

while true; do
  sleep 600
done
