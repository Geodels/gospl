import sys
import subprocess
import numpy as np
from scipy import spatial
from scipy import ndimage
from script import readOutput as output


# Specify most recent time in Ma BP
startMa = 0
# Specify deepest time in Ma BP
endMa = 20
# Specify paleo-displacements time interval in Ma
dtMa = 1

# Other parameters
res = 0.1
dt = 1.0e6
timeframe = np.arange(endMa, startMa, -dtMa)
factor = [
    0.2,
    0.4,
    0.6,
    0.8,
    1.0,
    0.1,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1.0,
    0.2,
    0.4,
    0.6,
    0.8,
    1.0,
]

backfile = "data/backward-elev0-20.npz"
outNb = 100  # Final time step for the model
sigma = 10  # Gaussian filter on tectonic grid

backMesh = np.load(backfile)

# Get initial input from Scotese
cmd0 = "python3 script/buildMesh.py -t=" + str(endMa) + " -d=data -s=100,30,15"
p = subprocess.Popen(cmd0, shell=True, stderr=subprocess.PIPE)
while True:
    out = p.stderr.read(1)
    if out == b"" and p.poll() is not None:
        break
    if out != "":
        sys.stdout.write(out.decode(sys.stdout.encoding))
        sys.stdout.flush()

# Main Loop...
for k in range(len(timeframe)):

    it = timeframe[k]

    ## Running gospl
    print("Running gospl")
    model1 = "inputs/model" + str(it) + "Ma.yml"
    cmd1 = "mpirun -np 4 python3 script/runModel.py -i " + model1
    p = subprocess.Popen(
        cmd1, shell=True, stdin=subprocess.DEVNULL, stderr=subprocess.PIPE
    )
    while True:
        out = p.stderr.read(1)
        if out == b"" and p.poll() is not None:
            break
        if out != "":
            sys.stdout.write(out.decode(sys.stdout.encoding))
            sys.stdout.flush()

    ## Scaled vertical tectonics from backward vs computed elevations
    fac = factor[k]
    print("Scale vertical tectoncis on regular grid")
    backstep = "z" + str(it - 1)
    modelinput = "inputs/model" + str(it) + "Ma.yml"
    tecfile = "input" + str(it) + "/vtec" + str(it) + "Ma"
    out = output.readOutput(filename=modelinput, step=outNb, uplift=False)
    out.buildLonLatMesh(res=res, nghb=3)
    if k == 0:
        shape = out.z.shape
        lon = np.linspace(0.0, 360, shape[1])
        lat = np.linspace(0.0, 180, shape[0])
        glon, glat = np.meshgrid(lon, lat)
        glonlat = np.dstack([glon.flatten(), glat.flatten()])[0]
        tree = spatial.cKDTree(glonlat, leafsize=10)
    print("Interpolate tectoncis on unstructured grid")
    tecto = fac * (backMesh[backstep] - out.z.flatten())
    tecto = ndimage.gaussian_filter(tecto, sigma)
    distances, indices = tree.query(out.lonlat, k=3)
    weights = np.divide(
        1.0, distances ** 2, out=np.zeros_like(distances), where=distances != 0
    )
    sumweights = np.sum(weights, axis=1)
    onIDs = np.where(sumweights == 0)[0]
    sumweights[sumweights == 0] = 1.0e-4
    tec = tecto.flatten() / dt
    uplift = np.sum(weights * tec[indices][:, :], axis=1) / sumweights
    if len(onIDs) > 0:
        uplift[onIDs] = tec[indices[onIDs, 0]]
    np.savez_compressed(tecfile, z=uplift)

    ## Running tectonically constrained gospl
    print("Running tectonically constrained gospl")
    model2 = "inputs/model" + str(it) + "Mat.yml"
    cmd2 = "mpirun -np 4 python3 script/runModel.py -i " + model2
    p = subprocess.Popen(
        cmd2, shell=True, stdin=subprocess.DEVNULL, stderr=subprocess.PIPE
    )
    while True:
        out = p.stderr.read(1)
        if out == b"" and p.poll() is not None:
            break
        if out != "":
            sys.stdout.write(out.decode(sys.stdout.encoding))
            sys.stdout.flush()

    ## Perform horizontal displacements and remeshing
    cmd3 = "python3 script/npzMesh.py -t=" + str(it - 1) + " -d=data -s=100,30,15 -i="
    if it > 16:
        cmd3 += model2 + " -n=100 -a=1 -r=20"
    elif it > 11:
        cmd3 += model2 + " -n=100 -a=1 -r=15"
    elif it > 6:
        cmd3 += model2 + " -n=100 -a=1 -r=10"
    else:
        cmd3 += model2 + " -n=100 -a=1 -r=5"

    p = subprocess.Popen(cmd3, shell=True, stderr=subprocess.PIPE)
    while True:
        out = p.stderr.read(1)
        if out == b"" and p.poll() is not None:
            break
        if out != "":
            sys.stdout.write(out.decode(sys.stdout.encoding))
            sys.stdout.flush()
