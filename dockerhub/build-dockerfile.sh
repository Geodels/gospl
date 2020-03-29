#!/bin/bash

# unset username
# unset password
# echo -n "username:"
# read username
# prompt="password:"
# while IFS= read -p "$prompt" -r -s -n 1 char
# do
#     if [[ $char == $'\0' ]]
#     then
#          break
#     fi
#     prompt='*'
#     password+="$char"
# done

# echo $'\n'

docker login 

echo -n "1. Do you wish to build the base gospl Docker container (y/n)?"
read answer
if [ "$answer" != "${answer#[Yy]}" ] ;then
    echo $'\n Build gospl base container... \n'
    docker build -t geodels/gospl-base:latest -f Dockerfile-debian .
    echo $'\n Push gospl base container... \n'
    docker push geodels/gospl-base:latest
fi

echo -n "2. Do you wish to build the main gospl Docker container (y/n)?"
read answer
if [ "$answer" != "${answer#[Yy]}" ] ;then
    echo $'\n Build gospl main container... \n'
    docker build -t geodels/gospl:latest -f Dockerfile .
    echo $'\n Push gospl main container... \n'
    docker push geodels/gospl:latest
fi

echo $'\n Docker images uploaded \n'
