## Global Scale LEM

> This repo is a specific branch of **gospl** tuned for HPC of global scale models.

#### Artemis HPC Sydney

```bash
module load gcc/4.9.3  python/3.6.5 petsc-gcc-mpich/3.11.1

export F77FLAGS=-fPIC
export FCFLAGS=-fPIC

python3 setup.py install --user
```

#### Performance tests

+ Resolution distribution for 10 million points (min:7.6 km, max:10.3 km, mean: 9.1 km)

| PTS | CPUS | WALLTIME | NODES |
| --- | --- | --- | --- |
| 10612062 | 8 | 17:08:58 | 1 |
| 10612062 | 16 | 08:58:32 | 2 |
| 10612062 | 32 | 04:57:58 | 4 |
| 10612062 | 64 | 03:32:15 | 6 |
| 10612062 | 96 | 02:41:17 | 6 |
| 10612062 | 128 | 02:40:57 | 7 |

+ Resolution distribution for 17 million points (min:4.8 km, max:7.6 km, mean: 6.0 km)

| PTS | CPUS | WALLTIME | NODES |
| --- | --- | --- | --- |
| 17004184 | 64 | 07:28:41 | 4 |
| 17004184 | 96 | 06:38:11 | 4 |
| 17004184 | 128 | 05:29:51 | 6 |
| 17004184 | 144 | 04:59:49 | 6 |
| 17004184 | 168 | 04:31:14 | 7 |
| 17004184 | 192 | 06:23:02 | 7 |
