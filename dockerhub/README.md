# Docker containers

Docker Container with all Dependencies

https://hub.docker.com/u/geodels/

A local script is given to perform the first 2 steps presented below (Docker base & Docker model).

You will run the following:
```
./build-dockerfile.sh
```

## Docker base

### Local installation

```
docker build -t geodels/gospl-base:latest -f Dockerfile-debian .
```

### Pushing the containers registry

```
docker login
```

```
docker push geodels/gospl-base:latest
```

## Docker model

### Local installation

```
docker build -t geodels/gospl:latest -f Dockerfile .
```

### Pushing the containers registry

```
docker login
```

```
docker push geodels/gospl:latest
```

## Docker paleoflow

### Local installation

```
docker build -t geodels/paleoflow:latest -f Dockerfile-paleoflow .
```

### Pushing the containers registry

```
docker login
```

```
docker push geodels/paleoflow:latest
```
