cd /live

export PATH=$PATH:/opt/mpich/bin/:/opt/petsc/bin
export SHELL=/bin/bash

jupyter-notebook --port=8888 --ip='0.0.0.0' --no-browser --allow-root --NotebookApp.token=''

while true; do
  sleep 600
done
