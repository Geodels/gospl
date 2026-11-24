import os
import numpy as np
import petsc4py
from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import getfacevelocity
    from gospl._fortran import advecupwind

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()


class WaveMesh(object):
    r"""
    Wave energy and wave shadows (opt-in ``waves:`` block,
    docs/DESIGN_CARBONATES.md 2.10, phase C7), on flat and global meshes.

    For each direction of the wave climate, a relative wave energy ``E``
    (0 to 1) travels across the sea along the local wave direction
    :math:`\hat{\mathbf{d}}`, grows over a fetch and is lost on land and in
    shallow water:

    .. math::
        \left(\hat{\mathbf{d}} \cdot \nabla + \frac{1}{L_f} + \kappa\right) E
        = \frac{1}{L_f}, \qquad
        \kappa = \frac{1}{L_d} + \frac{1}{L_s}
                  \max\left(0, 1 - \frac{w}{w_s}\right),

    with ``E = 0`` on land (identity rows: the shadow behind islands and
    headlands, from which the waves grow again over the fetch), ``w`` the
    water depth, ``L_f`` the fetch length over which the waves reach their
    full energy (``fetch``; without it there is no growth term), ``L_d`` the
    attenuation length in deep water (``decay``), ``w_s`` the depth below which
    the waves feel the bottom (``shoal_depth``) and ``L_s`` the attenuation
    length at the shoreline (``shoal_length``). On a flat mesh the domain
    edges facing the waves bring the open-ocean energy in (``E = 1``); a
    global mesh has no edges and needs the fetch.

    The local direction is a bearing (the direction the waves come from)
    turned into a tangent vector with the local east and north: one bearing
    everywhere, a map of bearings, or the ``zonal`` belts (trade winds from
    the east within 30 degrees of the equator, westerlies to 60 degrees,
    polar easterlies beyond). The climate is a list of directions with
    weights; ``E`` is their weighted sum. The operator is the first-order
    upwind finite-volume scheme of the orographic rain and of the tectonic
    advection (``getfacevelocity`` with a velocity per node +
    ``advecupwind``), whose coefficients are computed once per direction;
    the matrix is rebuilt every step because the coastline moves.

    The field is a response control (``wave``) of the carbonate factories and
    the energy of the coastal erosion (phase C4).
    """

    def __init__(self):
        """Cached coefficients per direction and the solver (registered in
        ``destroy_DMPlex``)."""
        self.waveE = None
        self._waveT = None
        self.coastEro = None
        self.coastBudget = None
        self._coastDet = None
        self._coastProv = None
        self._waveCoef = None
        self._waveMat = None
        self._ksp_wave = None
        if not getattr(self, "waveOn", False):
            return
        if not self.flatModel and self.waveFetch is None:
            raise ValueError(
                "waves: a global mesh has no domain edges to bring the waves in; "
                "set `fetch` (m), the length over which the waves grow to full "
                "energy over open water."
            )
        if self.flatModel and any(f == "zonal" for f, _ in self.waveClimate):
            raise ValueError(
                "waves: the `zonal` belts need latitudes, i.e. a global mesh; on a "
                "flat mesh give a bearing or a from_map.")
        self.waveE = np.zeros(self.lpoints, dtype=np.float64)
        if getattr(self, "coastOn", False):
            if not getattr(self, "seaDepo", False) or getattr(self, "nodep", False):
                raise ValueError(
                    "coastal_erosion: needs marine deposition (seadepo: True, "
                    "nodep: False): the eroded material is spread by the marine "
                    "routing.")
            self.coastEro = np.zeros(self.lpoints, dtype=np.float64)
            self._coastEroded = 0.0
            self._coastRubbleVol = 0.0
        return

    def _waveBearings(self, frm):
        """Bearing (degrees, the direction the waves come from) at every
        local node for one climate entry: a number, ``zonal`` belts, or a
        ``(file, key)`` map."""
        if isinstance(frm, tuple):
            d = np.load(frm[0] + ".npz")
            arr = np.asarray(d[frm[1]], dtype=np.float64)
            if arr.shape[0] != self.mpoints:
                raise ValueError("waves: bearing map %s has %d values, the mesh %d."
                                 % (list(frm), arr.shape[0], self.mpoints))
            return arr[self.locIDs]
        if frm == "zonal":
            lat = self._waveLatitude()
            b = np.where(np.abs(lat) < 30.0, 90.0, 270.0)       # trades, westerlies
            return np.where(np.abs(lat) >= 60.0, 90.0, b)       # polar easterlies
        return np.full(self.lpoints, float(frm))

    def _waveLatitude(self):
        """Latitude (degrees) of the local nodes (global mesh); on a flat mesh
        the ``y`` coordinate is not a latitude, so the zonal belts are refused
        there by the parser."""
        xyz = self.lcoords
        r = np.linalg.norm(xyz, axis=1)
        return np.degrees(np.arcsin(np.clip(xyz[:, 2] / np.maximum(r, 1.0), -1.0, 1.0)))

    def _waveVelocity(self, bearing):
        """Unit wave-travel vectors per node from bearings (the direction the
        waves come from, 0 = north, 90 = east): along ``(x, y)`` on a flat
        mesh, in the local tangent plane (east, north) on a sphere."""
        a = np.deg2rad(bearing)
        ue, un = -np.sin(a), -np.cos(a)                     # travel away from `from`
        vel = np.zeros((self.lpoints, 3))
        if self.flatModel:
            vel[:, 0] = ue
            vel[:, 1] = un
            return vel
        xyz = self.lcoords
        lon = np.arctan2(xyz[:, 1], xyz[:, 0])
        lat = np.arcsin(np.clip(xyz[:, 2] / np.maximum(np.linalg.norm(xyz, axis=1), 1.0),
                                -1.0, 1.0))
        east = np.column_stack([-np.sin(lon), np.cos(lon), np.zeros_like(lon)])
        north = np.column_stack([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon),
                                 np.cos(lat)])
        return ue[:, None] * east + un[:, None] * north

    def _waveInflow(self, vel, grazing=0.05):
        """Domain-edge nodes (flat mesh) where the waves enter: the local
        travel direction points into the domain, ``dot(v, n) < -grazing`` with
        ``n`` the outward edge normal estimated from the neighbours (the test
        of ``_oroInflowBorders``, with a direction per node). Rank-local."""
        b = np.asarray(self.advectBorders, dtype=np.int64)
        if b.size == 0:
            return b
        ngb = self.FVmesh_ngbID[b]
        valid = ngb >= 0
        xy = self.lcoords[:, :2]
        cnt = np.maximum(valid.sum(axis=1), 1)
        cx = np.where(valid, xy[np.maximum(ngb, 0), 0], 0.0).sum(axis=1) / cnt
        cy = np.where(valid, xy[np.maximum(ngb, 0), 1], 0.0).sum(axis=1) / cnt
        nx, ny = xy[b, 0] - cx, xy[b, 1] - cy
        nrm = np.hypot(nx, ny)
        cosv = np.divide(vel[b, 0] * nx + vel[b, 1] * ny, nrm, out=np.zeros_like(nrm),
                         where=nrm > 0.0)
        return b[cosv < -grazing]

    def _waveBasis(self):
        """Local tangent basis ``(e1, e2)`` per node: ``x, y`` on a flat mesh,
        east and north on a sphere."""
        if self.flatModel:
            e1 = np.zeros((self.lpoints, 3))
            e2 = np.zeros((self.lpoints, 3))
            e1[:, 0] = 1.0
            e2[:, 1] = 1.0
            return e1, e2
        xyz = self.lcoords
        lon = np.arctan2(xyz[:, 1], xyz[:, 0])
        lat = np.arcsin(np.clip(xyz[:, 2] / np.maximum(np.linalg.norm(xyz, axis=1), 1.0),
                                -1.0, 1.0))
        e1 = np.column_stack([-np.sin(lon), np.cos(lon), np.zeros_like(lon)])
        e2 = np.column_stack([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon),
                              np.cos(lat)])
        return e1, e2

    def _waveGradient(self, f):
        """
        Least-squares gradient of ``f`` in the local tangent plane, as a 3-D
        vector per node: the plane fit of each node's neighbours (exact for a
        linear field). Owned rows are complete; ghost rows are taken from
        their owners (collective). Scratch: ``tmp``, ``tmpL``.
        """
        e1, e2 = self._waveBasis()
        ngb = self.FVmesh_ngbID
        ok = ngb >= 0
        nb = np.where(ok, ngb, np.arange(self.lpoints)[:, None])
        dx = self.lcoords[nb] - self.lcoords[:, None, :]
        a = np.einsum("nkc,nc->nk", dx, e1) * ok
        b = np.einsum("nkc,nc->nk", dx, e2) * ok
        df = (f[nb] - f[:, None]) * ok
        saa, sab, sbb = (a * a).sum(1), (a * b).sum(1), (b * b).sum(1)
        saf, sbf = (a * df).sum(1), (b * df).sum(1)
        det = saa * sbb - sab * sab
        good = np.abs(det) > 1.0e-30
        ga = np.divide(sbb * saf - sab * sbf, det, out=np.zeros(self.lpoints), where=good)
        gb = np.divide(saa * sbf - sab * saf, det, out=np.zeros(self.lpoints), where=good)
        for arr in (ga, gb):
            self.tmpL.setArray(arr)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.dm.globalToLocal(self.tmp, self.tmpL)
            arr[:] = self.tmpL.getArray()
        return ga[:, None] * e1 + gb[:, None] * e2

    def _waveRefract(self, vel0, depth):
        """
        Refracted wave directions (refraction, phase C7b): deeper than
        ``shoal_depth`` the offshore direction ``vel0``; shallower, turned
        towards the normal of the depth contours (pointing to shallower water)
        so that ``sin(theta) = sqrt(w / w_s) sin(theta_0)``, Snell's law with
        the shallow-water celerity ``sqrt(g w)`` (``theta`` measured from the
        contour normal). Waves travelling away from the shore are not turned.
        The depth may be smoothed over ``refraction_length`` first.
        Collective (gradient sync, smoothing solve).
        """
        w = np.maximum(depth, 0.0)
        if self.waveRefLength is not None:
            w = self._waveSmooth(w)
        g = self._waveGradient(w)
        gn = np.linalg.norm(g, axis=1)
        n = np.divide(-g, gn[:, None], out=np.zeros_like(g), where=gn[:, None] > 0.0)
        cos0 = np.einsum("nc,nc->n", vel0, n)
        t = vel0 - cos0[:, None] * n
        tn = np.linalg.norm(t, axis=1)
        t = np.divide(t, tn[:, None], out=np.zeros_like(t), where=tn[:, None] > 0.0)
        sin0 = tn
        ratio = np.sqrt(np.clip(w / self.waveShoalDepth, 0.0, 1.0))
        sin1 = np.clip(ratio * sin0, 0.0, 1.0)
        cos1 = np.sqrt(1.0 - sin1 * sin1)
        vel = cos1[:, None] * n + sin1[:, None] * t
        turn = (cos0 > 0.0) & (gn > 0.0) & (depth > 0.0) & (depth < self.waveShoalDepth)
        return np.where(turn[:, None], vel, vel0)

    def _waveSmooth(self, f):
        """``f`` smoothed over ``refraction_length``, (I - L^2 lap) s = f
        (the open-water index operator, built for this length). Collective."""
        if getattr(self, "_waveSmMat", None) is None:
            kd = np.full(self.lpoints, self.waveRefLength ** 2, dtype=np.float64)
            self._waveSmMat = self._buildDiffMat(kd)
            ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
            ksp.setType("fgmres")
            ksp.getPC().setType("bjacobi")
            ksp.setTolerances(rtol=1.0e-10, max_it=2000)
            ksp.setOptionsPrefix("wavesm_")
            petsc4py.PETSc.Options()["wavesm_sub_pc_factor_shift_type"] = "nonzero"
            ksp.setFromOptions()
            ksp.setOperators(self._waveSmMat, self._waveSmMat)
            self._ksp_wavesm = ksp
        self.tmpL.setArray(f)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.tmp1.set(0.0)
        self._ksp_wavesm.solve(self.tmp, self.tmp1)
        self.dm.globalToLocal(self.tmp1, self.tmpL)
        return self.tmpL.getArray().copy()

    def _waveCoefficients(self):
        """Upwind advection coefficients (unit speed, so the operator is
        ``dot(d, grad)`` in 1/m) and the inflow edges (flat meshes) of every
        direction of the climate. Computed once (the mesh does not change);
        ``getfacevelocity`` sets module state of the Fortran kernels, so the
        coefficients are copied at once."""
        self._waveCoef = []
        depth = self.sealevel - self.hLocal.getArray()
        if getattr(self, "_waveBear", None) is None:
            self._waveBear = [self._waveBearings(frm) for frm, _ in self.waveClimate]
        self._waveVel = []
        for k, (frm, _) in enumerate(self.waveClimate):
            bearing = self._waveBear[k]
            vel = self._waveVelocity(bearing)
            if self.waveRefraction:
                vel = self._waveRefract(vel, depth)
            self._waveVel.append(vel)
            getfacevelocity(self.lpoints, vel)
            lcoeff = advecupwind(self.lpoints, 1.0)
            inflow = self._waveInflow(vel) if self.flatModel else np.zeros(0, dtype=np.int64)
            self._waveCoef.append((lcoeff[:, 0] - 1.0, lcoeff[:, 1:].copy(), inflow))

    def _waveSolve(self, diag, off, fixed, rhs):
        """Assemble ``diag`` + ``off`` with identity rows at ``fixed`` and
        solve with the right-hand side ``rhs`` (the imposed values on the
        fixed rows, the source elsewhere); returns the local solution.
        Collective."""
        d = diag.copy()
        o = off.copy()
        d[fixed] = 1.0
        o[fixed, :] = 0.0
        coeffs = np.empty((self.lpoints, 1 + self.maxnb), dtype=np.float64)
        coeffs[:, 0] = d
        coeffs[:, 1:] = o[:, : self.maxnb]
        if self._waveMat is not None:
            self._waveMat.destroy()
        self._waveMat = self._assembleDiffMatCSR(coeffs)
        if self._ksp_wave is None:
            ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
            ksp.setType("gmres")
            ksp.getPC().setType("bjacobi")
            ksp.setTolerances(rtol=1.0e-10, atol=1.0e-14, max_it=2000)
            ksp.setOptionsPrefix("wave_")
            petsc4py.PETSc.Options()["wave_sub_pc_factor_shift_type"] = "nonzero"
            ksp.setFromOptions()
            self._ksp_wave = ksp
        self._ksp_wave.setOperators(self._waveMat, self._waveMat)
        self.tmpL.setArray(rhs)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.tmp1.set(0.0)
        self._ksp_wave.solve(self.tmp, self.tmp1)
        if self._ksp_wave.getConvergedReason() < 0 and MPIrank == 0:
            print("[waves] wave energy solve did not converge (%d)"
                  % self._ksp_wave.getConvergedReason(), flush=True)
        self.dm.globalToLocal(self.tmp1, self.tmpL)
        return self.tmpL.getArray().copy()

    def waveEnergy(self):
        """
        Relative wave energy ``waveE`` (0 to 1) of the current surface: the
        weighted sum over the wave climate of the solutions of
        ``(dot(d, grad) + 1/L_f + kappa) E = 1/L_f`` (see the class). Collective. Scratch:
        ``tmp``, ``tmpL``, ``tmp1``.
        """
        if self._waveT is not None and self._waveT == self.tNow and self.waveE is not None:
            return self.waveE
        t0 = process_time()
        if self._waveCoef is None or self.waveRefraction:
            self._waveCoefficients()
        h = self.hLocal.getArray()
        w = self.sealevel - h
        land = w <= 0.0
        kappa = 1.0 / self.waveDecay + np.maximum(0.0, 1.0 - w / self.waveShoalDepth) \
            / self.waveShoalLength
        grow = 0.0 if self.waveFetch is None else 1.0 / self.waveFetch
        E = np.zeros(self.lpoints, dtype=np.float64)
        for (frm, wgt), (adiag, off, inflow) in zip(self.waveClimate, self._waveCoef):
            fixed = land.copy()
            # the growth source over open water, and the open-ocean energy on
            # the inflow edges (flat meshes)
            val = np.where(land, 0.0, grow)
            inn = np.zeros(self.lpoints, dtype=bool)
            inn[inflow] = True
            inn &= ~land
            fixed |= inn
            val[inn] = 1.0
            Ek = self._waveSolve(adiag + kappa + grow, off, fixed, val)
            Ek[land] = 0.0
            E += wgt * np.clip(Ek, 0.0, 1.0)
        self.waveE = E
        self._waveT = self.tNow
        if MPIrank == 0 and self.verbose:
            print("[waves] wave energy (%0.02f seconds)" % (process_time() - t0),
                  flush=True)
        return E

    def _coastToRubble(self, litv, det, fine):
        """Route the eroded reef framework as the ``rubble`` lithology: moves
        column ``reef`` of ``litv`` (deposit-equivalent volumes per lithology,
        m^3, modified in place) to the rubble column at the same solid volume,
        ``V_r (1 - phi0_reef) / (1 - phi0_rubble)``, and corrects the routed
        bulk ``det`` and fine-group ``fine`` volumes. Returns ``(det, fine,
        solid)`` with ``solid`` the reef solid volume converted per node (m^3).
        Rank-local."""
        ir, it = self.coastRubble
        vr = litv[:, ir].copy()
        solid = vr * (1.0 - self.litPhi0[ir])
        conv = solid / (1.0 - self.litPhi0[it])
        litv[:, ir] = 0.0
        litv[:, it] += conv
        det = det + conv - vr
        if self.stratLith:
            fine = fine + (conv if self.litFine[it] else 0.0) \
                - (vr if self.litFine[ir] else 0.0)
        return det, fine, solid

    def coastalErosion(self):
        """
        Wave erosion of the coast and the shallow seafloor (phase C4), after
        Pastier et al. (2019): with the wave energy ``E`` (``waveEnergy``) and
        the erodibility ``K`` of the exposed layer (stratigraphic ``stratK`` of
        the top layer times the lithology and duricrust factor),

        - a marine cell at depth ``w`` erodes ``rate * E * exp(-w /
          depth_scale) * K * dt``;
        - a land cell up to ``cliff_height`` above sea level with a sea
          neighbour (a cliff) erodes ``rate * E_sea * K * dt`` (``E_sea`` the
          largest energy of its sea neighbours), at most down to sea level: the
          cliff retreats and leaves a platform, abraded in turn;

        both capped at ``max_erosion`` per step. The eroded material is taken
        from the stratigraphy (``erodeStrat``, its deposit-equivalent volume and
        composition) and handed to this step's marine routing (``seaChange``),
        like the loose carbonate. With ``rubble`` the eroded reef framework is
        routed as that lithology (reef rubble, usually ``carb_sand``): the same
        solid volume, re-expressed at the rubble's surface porosity, so the
        routed deposit-equivalent volume changes by the porosity ratio.
        Elevation, cumulative erosion and the rate field are updated. Collective (wave solve, halo sync, budget).
        Scratch: ``tmp``, ``tmpL``.
        """
        E = self.waveEnergy()
        h = self.hLocal.getArray().copy()
        w = self.sealevel - h
        kf = self._surfaceK() * self._surfaceLithoK() if self.stratNb > 0 \
            else np.ones(self.lpoints)
        kf = kf * np.ones(self.lpoints)
        sea = w > 0.0
        ero = np.where(sea, self.coastRate * E * np.exp(-np.maximum(w, 0.0) / self.coastDepth),
                       0.0)
        # cliffs: land cells next to the sea, low enough above it
        ngb = self.FVmesh_ngbID
        ok = ngb >= 0
        nb = np.where(ok, ngb, np.arange(self.lpoints)[:, None])
        nsea = ok & sea[nb]
        esea = np.where(nsea, E[nb], 0.0).max(axis=1)
        cliff = (~sea) & nsea.any(axis=1) & (h - self.sealevel <= self.coastCliff)
        ero = np.where(cliff, self.coastRate * esea, ero)
        ero = ero * kf * self.dt
        # a cliff is cut down to sea level at most
        ero = np.where(cliff, np.minimum(ero, np.maximum(h - self.sealevel, 0.0) + 1.0e-6),
                       ero)
        ero = np.clip(ero, 0.0, self.coastMax)
        # halo: the owners' values on the ghosts
        self.tmpL.setArray(ero)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        ero = self.tmpL.getArray().copy()

        # remove it from the stratigraphy (keep the river step's arrays)
        lit = getattr(self, "_litMulti", False)
        vol = ero * self.larea
        if self.stratNb > 0:
            saved = {k: getattr(self, k).copy() for k in
                     ("thCoarse", "thFine", "provEro", "thL")
                     if getattr(self, k, None) is not None}
            self.tmpL.setArray(-ero)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.erodeStrat()
            det = self.thCoarse * self.dt * self.larea
            fine = np.zeros(self.lpoints)
            if self.stratLith:
                fine = self.thFine * self.dt * self.larea
                det = det + fine
            litv = None
            if lit:
                litv = np.column_stack([self.thL[l] * self.dt * self.larea
                                        for l in range(self.litNb)])
            # provenance of what was eroded (per-class volumes, sum == det)
            provc = None
            if getattr(self, "provOn", False) and getattr(self, "provEro", None) is not None:
                provc = self.provEro * self.dt * self.larea[:, None]
            for k, v in saved.items():
                setattr(self, k, v)
        else:
            det, fine, litv, provc = vol.copy(), np.zeros(self.lpoints), None, None
        rubble = 0.0
        if litv is not None and getattr(self, "coastRubble", None) is not None:
            det0 = det.copy()
            det, fine, rubble = self._coastToRubble(litv, det, fine)
            if provc is not None:
                # same composition, volume re-expressed at the rubble porosity
                provc = provc * np.divide(det, det0, out=np.ones_like(det),
                                          where=det0 > 0.0)[:, None]
        self._coastDet = det
        self._coastFine = fine
        self._coastL = litv
        self._coastProv = provc

        # elevation, cumulative erosion, rate
        self.tmpL.setArray(-ero)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.cumED.axpy(1.0, self.tmp)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)
        self.tmpL.setArray(-ero / self.dt)
        self.EbLocal.axpy(1.0, self.tmpL)
        self.coastEro += ero

        own = self.inIDs == 1
        rub = rubble[own].sum() if np.ndim(rubble) else 0.0
        b = MPI.COMM_WORLD.allreduce(np.array([vol[own].sum(), det[own].sum(), rub]),
                                     op=MPI.SUM)
        # (bulk eroded, routed deposit-equivalent, reef solid routed as rubble)
        self.coastBudget = b
        self._coastEroded += b[0]
        self._coastRubbleVol += b[2]
        if MPIrank == 0 and self.verbose:
            print("[coast] eroded %0.4g km3 (routed %0.4g km3)" % (b[0] * 1e-9, b[1] * 1e-9),
                  flush=True)
        return
