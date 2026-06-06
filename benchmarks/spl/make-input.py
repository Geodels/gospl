import os
import numpy as np
import meshio
import xarray as xr
from scripts import umeshFcts as ufcts

# =============================================================================
# SPL Benchmark — Mesh Generation
# =============================================================================
# Creates a planar unstructured mesh for a single-basin Stream Power Law
# benchmark with uniform uplift and a single central outlet.
#
# Domain   : nx × ny cells at dx resolution  →  20 km × 20 km
# Outlet   : centre of the southern edge (fixed at z = 0)
# Uplift   : uniform U everywhere except the southern boundary row
# Output   : gospl_mesh.npz  (vertex, cells, elevation, uplift)
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Domain parameters
# -----------------------------------------------------------------------------
nx, ny = 200, 200       # number of nodes in x and y
dx     = 100.0          # node spacing (m)

x = np.arange(nx) * dx  # x coordinates (m)  0 → 19900 m
y = np.arange(ny) * dx  # y coordinates (m)  0 → 19900 m
X, Y = np.meshgrid(x, y)

print(f"Domain        : {nx*dx/1e3:.0f} km × {ny*dx/1e3:.0f} km")
print(f"Resolution    : {dx:.0f} m")
print(f"Nodes         : {nx*ny:,}")

# -----------------------------------------------------------------------------
# 2. Initial topography
# -----------------------------------------------------------------------------
# Elevation increases linearly away from the outlet (gentle radial ramp)
# with small random perturbations to seed drainage network organisation.
# The outlet node is pinned to z = 0 (base level).

outlet_i = nx // 2      # outlet column — centre of southern edge
outlet_j = ny - 1       # outlet row    — southern boundary

# Distance from each node to the outlet (in node units)
r = np.sqrt((X/dx - outlet_i)**2 + (Y/dx - outlet_j)**2)

np.random.seed(42)                          # reproducible perturbations
z0  = 0.01 * r                              # gentle radial slope (m)
z0 += 0.1 * np.random.randn(ny, nx)        # small noise to seed channels
z0[outlet_j, outlet_i] = 0.0               # pin outlet to base level

print(f"\nInitial topography:")
print(f"  Elevation range : {z0.min():.2f} – {z0.max():.2f} m")
print(f"  Outlet node     : ({outlet_i}, {outlet_j})  z = 0")

# -----------------------------------------------------------------------------
# 3. Uplift field
# -----------------------------------------------------------------------------
# Uniform uplift U everywhere except the southern boundary row,
# which is held fixed at z = 0 (open boundary / base level).

U      = 4.e-4                      # uplift rate (m/yr) = 0.4 mm/yr
uplift = np.full_like(z0, U)
uplift[-1, :] = 0.0                 # no uplift on southern boundary

print(f"\nUplift field:")
print(f"  U (interior)    : {U:.1e} m/yr  ({U*1e3:.2f} mm/yr)")
print(f"  U (S boundary)  : 0.0  (fixed base level)")

# -----------------------------------------------------------------------------
# 4. Build xarray dataset for mesh generation
# -----------------------------------------------------------------------------
ds = xr.Dataset({
    'elev': xr.DataArray(
        data   = z0,
        dims   = ['y', 'x'],
        coords = {'x': x, 'y': y},
        attrs  = {'units': 'm', 'long_name': 'initial surface elevation'}
    ),
    'tec': xr.DataArray(
        data   = uplift,
        dims   = ['y', 'x'],
        coords = {'x': x, 'y': y},
        attrs  = {'units': 'm/yr', 'long_name': 'uplift rate'}
    ),
})
ds['cellwidth'] = (['y', 'x'], dx * np.ones((ny, nx)))

# -----------------------------------------------------------------------------
# 5. Generate planar unstructured mesh
# -----------------------------------------------------------------------------
output_path = "meshing"

if os.path.exists(output_path):
    os.system(f"rm -rf {output_path}")
os.makedirs(output_path)

print(f"\nBuilding unstructured mesh → {output_path}/")
ufcts.planarMesh(ds, output_path, fvtk='planar.vtk', fumpas=True, voro=True)

# -----------------------------------------------------------------------------
# 6. Interpolate fields onto unstructured mesh
# -----------------------------------------------------------------------------
ufile    = os.path.join(output_path, 'base2D.nc')
var_name = 'data'
mapds    = xr.open_dataset(ufile)

print("Interpolating elevation and uplift onto unstructured mesh nodes...")
ufcts.inter2UGRID(ds[['elev', 'tec']], mapds, output_path,
                  var_name, type='face', latlon=False)
data_ds = xr.open_dataset(os.path.join(output_path, var_name + '.nc'))

# -----------------------------------------------------------------------------
# 7. Mesh diagnostics
# -----------------------------------------------------------------------------
n_nodes  = mapds.sizes['nCells']
ucoords  = np.column_stack([
    mapds['xCell'].values,
    mapds['yCell'].values,
    mapds['zCell'].values,
])
ufaces   = mapds['cellsOnVertex'].values - 1  # (nVertices, 3) — 0-based

dcEdge    = mapds['dcEdge'].values
edge_min  = np.round(dcEdge.min()  / 1e3, 2)
edge_max  = np.round(dcEdge.max()  / 1e3, 2)
edge_mean = np.round(dcEdge.mean() / 1e3, 2)

print(f"\nUnstructured mesh:")
print(f"  Nodes           : {n_nodes:,}")
print(f"  Triangles       : {ufaces.shape[0]:,}")
print(f"  Edge length     : min {edge_min} km | "
      f"mean {edge_mean} km | max {edge_max} km")

# -----------------------------------------------------------------------------
# 8. Read triangulation and save goSPL input file
# -----------------------------------------------------------------------------
mesh   = meshio.read(os.path.join(output_path, 'planar.vtk'))
vertex = mesh.points
cells  = mesh.cells_dict['triangle']

print(f"\ngoSPL mesh:")
print(f"  Vertices        : {len(vertex):,}")
print(f"  Triangles       : {len(cells):,}")
print(f"  Elevation range : {data_ds.elev.data.min():.2f} – "
      f"{data_ds.elev.data.max():.2f} m")

meshname = "gospl_mesh"
np.savez_compressed(
    meshname,
    v = vertex,         # (N, 3)  node coordinates
    c = cells,          # (Nc, 3) triangle connectivity  (1-based)
    z = data_ds.elev.data,   # (N,)    initial elevation (m)
    t = data_ds.tec.data,    # (N,)    uplift rate (m/yr)
)
print(f"\nSaved → {meshname}.npz")
print(f"  Arrays: v (vertices), c (cells), z (elevation), t (uplift)")

# -----------------------------------------------------------------------------
# 9. Cleanup
# -----------------------------------------------------------------------------
os.system(f"rm -rf {output_path}")
print(f"\nTemporary meshing directory removed.")
print(f"\nReady to run goSPL.  Point your YAML to: {meshname}.npz")