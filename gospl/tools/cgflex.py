import os
import gc
import sys
import glob
import pickle
import pyproj
import petsc4py
import numpy as np

import xarray as xr
from mpi4py import MPI
from pathlib import Path
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gflex.f2d import F2D

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD


class GlobalFlex(object):
    """
    This class defines global flexural responses.

    """

    def __init__(self):
        """
        The initialisation of `globalFlex` class.
        """

        if self.gflexOn:
            
            self.cumEDFlex = self.hLocal.duplicate()
            self.cumEDFlex.set(0.0)
            self.cumFlexL = self.hLocal.duplicate()
            self.cumFlexL.set(0.0)
            interpData = np.load(self.Interp)
            self.gflex_ids = interpData["ids"]
            self.gflex_oIDs = interpData["oids"]
            self.gflex_wghts = interpData["wghts"]
            self.gflex_sumwght = interpData["sumwght"]
            self.nextFlex = self.tNow + self.gflexStep

            if MPIrank == 0:        
                interpRData = np.load(self.rInterp)
                self.rflex_wghts = interpRData["wghts"]
                self.rflex_ids = interpRData["ids"]
                self.rflex_oIDs = interpRData["oids"]
                self.rflex_denum = interpRData["denum"]
                del interpRData
                gc.collect()

            # All this is hard coded to produce a 0.25 deg resolution grid globally
            self.res = 0.25
            self.regshape = (721, 1441)
            self.long = np.linspace(-180.0, 180.0, self.regshape[1])
            self.lati = np.linspace(-90.0, 90.0, self.regshape[0])
            del interpData
            gc.collect()
        else:
            self.nextFlex = self.tEnd + 1.e6

        return
    
    def _makeBnds(self, dlon=20, dlat=15):

        lonrange = np.arange(-140,140+dlon,dlon)
        latrange = np.arange(-75,75+dlat,dlat)

        listbnds = []
        for lon in range(len(lonrange)-1):
            min_lon = lonrange[lon]
            max_lon = lonrange[lon+1]
            for lat in range(len(latrange)-1):
                min_lat = latrange[lat]
                max_lat = latrange[lat+1]
                listbnds.append([min_lon,max_lon,min_lat,max_lat])
        
        del lonrange, latrange
        gc.collect()
        
        return listbnds

    def _shiftValues(self, ds):

        zval = ds.elevation.values.copy()
        dval = ds.erodep.values.copy()
        lx = int((zval.shape[1]-1)/2)
        shifth = np.zeros(zval.shape) 
        shifth[:,:lx+1] = zval[:,lx:]
        shifth[:,lx+1:] = zval[:,:lx]
        shifted = np.zeros(dval.shape) 
        shifted[:,:lx+1] = dval[:,lx:]
        shifted[:,lx+1:] = dval[:,:lx]

        ds['shiftz'] = (('latitude','longitude'),shifth)
        ds['shiftz'].attrs['units'] = 'metres'

        ds['shifted'] = (('latitude','longitude'),shifted)
        ds['shifted'].attrs['units'] = 'metres'

        del zval, dval, shifth, shifted
        gc.collect()

        return ds

    def _interpUTM(self, ds):

        xmax,xmin = ds.x.values.max(),ds.x.values.min()
        ymax,ymin = ds.y.values.max(),ds.y.values.min()
        dx = (xmax-xmin)/len(ds.x.values)
        dy = (ymax-ymin)/len(ds.y.values)
        rdx = np.round(dx, -3)
        rdy = np.round(dy, -3)

        rxmin = np.round(xmin, -3)
        if rxmin < xmin:
            rxmin += 1000
        rxmax = np.round(xmax, -3)
        if rxmax > xmax:
            rxmax -= 1000

        rymin = np.round(ymin, -3)
        if rymin < ymin:
            rymin += 1000
        rymax = np.round(ymax, -3)
        if rymax > ymax:
            rymax -= 1000
        newx = np.arange(rxmin,xmax,rdx)
        newy = np.arange(rymin,ymax,rdy)

        return ds.interp(x=newx, y=newy, method="cubic")

    def _makeUTM(self, ds, lonlat_bnds, dll=10):

        min_lon = lonlat_bnds[0]-5
        max_lon = lonlat_bnds[1]+5
        min_lat = lonlat_bnds[2]-5
        max_lat = lonlat_bnds[3]+5
        ds.rio.write_crs(4326, inplace=True)
        
        lon_bnds = [min_lon, max_lon]
        lat_bnds = [min_lat, max_lat]
        ds_clip = ds.sel(latitude=slice(*lat_bnds), longitude=slice(*lon_bnds))
        
        ds_srtm = ds_clip.rename({'longitude': 'x','latitude': 'y'})
        epsg = ds_srtm.rio.estimate_utm_crs()
        clip_utm = ds_srtm.rio.reproject(epsg)
        clip_utm = clip_utm.sortby(clip_utm.x)
        clip_utm = clip_utm.sortby(clip_utm.y)
        x_bnds = [clip_utm.x.min(), clip_utm.x.max()]
        y_bnds = [clip_utm.y.min(), clip_utm.y.max()]

        lon_bnds2 = [min_lon-dll, max_lon+dll]
        lat_bnds2 = [min_lat-dll, max_lat+dll]
        ds_clip2 = ds.sel(latitude=slice(*lat_bnds2), longitude=slice(*lon_bnds2))
        clip_utm2 = ds_clip2.rio.reproject(epsg)
        clip_utm2 = clip_utm2.sortby(clip_utm2.x)
        clip_utm2 = clip_utm2.sortby(clip_utm2.y)

        ds2_clip = clip_utm2.sel(x=slice(*x_bnds), y=slice(*y_bnds))
        ds2_clip = ds2_clip.where(ds2_clip.elevation<20000)

        interp = self._interpUTM(ds2_clip)
        
        ds_srtm.close()
        ds_clip.close()
        ds2_clip.close()
        ds_clip2.close()
        clip_utm2.close()
        gc.collect()

        return interp 

    def _makeLonLat(self, utmds, bnds, res):

        test = utmds.copy()
        test.rio.write_crs(test.rio.crs, inplace=True)
        testll = test.rio.reproject('EPSG:4326')
        testll = testll.sortby(testll.x)
        testll = testll.sortby(testll.y)
        testll = testll.rename({'x': 'longitude','y': 'latitude'})

        lon_bnds = [bnds[0]-res*2, bnds[1]+res*2]
        lat_bnds = [bnds[2]-res*2, bnds[3]+res*2]
        map = testll.sel(latitude=slice(*lat_bnds), longitude=slice(*lon_bnds))
        map = map.where(map.flex<200000)
        
        test.close()

        return map

    def _interpLonLat(self, ds, bnds, res):
        rxmin,rxmax = bnds[0], bnds[1]
        rymin,rymax = bnds[2], bnds[3]
        newx = np.arange(rxmin,rxmax+res,res)
        newy = np.arange(rymin,rymax+res,res)
        return ds.interp(longitude=newx, latitude=newy, method="nearest")

    def _interpLonLat2(self, ds, bnds, res):
        rxmin,rxmax = bnds[0], bnds[1]
        rymin,rymax = bnds[2], bnds[3]
        newx = np.arange(rxmin+res,rxmax,res)
        newy = np.arange(rymin,rymax+res,res)
        return ds.interp(longitude=newx, latitude=newy, method="nearest")

    def _back2LonLat(self, utm_ds, bnds, res):
        tmp = self._makeLonLat(utm_ds,bnds,res)
        return self._interpLonLat(tmp,bnds,res)

    def _checkData(self, ds, bnds):
        lon_bnds = [bnds[0], bnds[1]]
        lat_bnds = [bnds[2], bnds[3]]
        return ds.sel(latitude=slice(*lat_bnds), longitude=slice(*lon_bnds))

    def _getFlex(self, ds, fds, fdse):

        nfds = xr.Dataset({
            'empty': xr.DataArray(np.nan, coords=dict(latitude=ds.latitude.values, 
                                            longitude=ds.longitude.values), 
                            dims=("latitude", "longitude")) 
            }
            )

        shape = nfds.empty.values.shape
        sflexe = fdse.flex.values.copy()

        newflex = np.zeros(shape)
        newflex[:,1:-1] = fds.flex.values.copy()
        
        datalat_min,datalat_max = fdse.latitude.values.min(),fdse.latitude.values.max()

        ilat_min = list(nfds.latitude.values).index(nfds.sel(latitude=datalat_min, method='nearest').latitude)
        ilat_max = list(nfds.latitude.values).index(nfds.sel(latitude=datalat_max, method='nearest').latitude)
        ilon140 = list(nfds.longitude.values).index(nfds.sel(longitude=140.0, method='nearest').longitude)
        nlon140 = list(nfds.longitude.values).index(nfds.sel(longitude=-140.0, method='nearest').longitude)

        ilon0 = list(fdse.longitude.values).index(fdse.sel(longitude=0.0, method='nearest').longitude)
        ilon40 = list(fdse.longitude.values).index(fdse.sel(longitude=40.0, method='nearest').longitude)
        nlon40 = list(fdse.longitude.values).index(fdse.sel(longitude=-40.0, method='nearest').longitude)
        
        newflex[ilat_min:ilat_max+1,0:nlon140+1] =  sflexe[:,ilon0:ilon40+1]
        newflex[ilat_min:ilat_max+1,ilon140:] =  sflexe[:,nlon40:ilon0+1]
        newflex[0:ilat_min,0] = newflex[0:ilat_min,1]
        newflex[0:ilat_min,-1] = newflex[0:ilat_min,-2]
        newflex[ilat_max+1:,0] = newflex[ilat_max+1:,1]
        newflex[ilat_max+1:,-1] = newflex[ilat_max+1:,-2]

        nfds['flex'] = (('latitude','longitude'),newflex)
        nfds['flex'].attrs['units'] = 'metres'

        return nfds[['flex']]

    def _lon_lat_box(self, lon_bounds, lat_bounds, refinement=2):    
        lons = []
        lats = []
        lons.append(np.linspace(lon_bounds[0], lon_bounds[-1], num=refinement))
        lats.append(np.linspace(lat_bounds[0], lat_bounds[0], num=refinement))
        lons.append(np.linspace(lon_bounds[-1], lon_bounds[-1], num=refinement))
        lats.append(np.linspace(lat_bounds[0], lat_bounds[-1], num=refinement))                
        lons.append(np.linspace(lon_bounds[-1], lon_bounds[0], num=refinement))
        lats.append(np.linspace(lat_bounds[-1], lat_bounds[-1], num=refinement))
        lons.append(np.linspace(lon_bounds[0], lon_bounds[0], num=refinement))
        lats.append(np.linspace(lat_bounds[-1], lat_bounds[0], num=refinement))
        return np.concatenate(lons), np.concatenate(lats)

    def _getSouthPole(self, ds):

        # Polar Stereographic South (71S,0E) epsg:3031
        source_crs = 'epsg:3031' 
        target_crs = 'epsg:4326'
        latlon_to_spolar = pyproj.Transformer.from_crs(target_crs, source_crs)
        sblon, sblat = self._lon_lat_box([-180,180],[-90,-60], refinement=100)
        sbx, sby = latlon_to_spolar.transform(sblat, sblon)

        sxmin,sxmax = sbx.min(),sbx.max()
        symin,symax = sby.min(),sby.max()

        southpole = ds.sel(latitude=slice(*[-90, -35]))
        southpole.rio.write_crs(target_crs, inplace=True)
        sPolar = southpole.rio.reproject(source_crs)
        sPolar = sPolar.sortby(sPolar.x)
        sPolar = sPolar.sortby(sPolar.y)
        southClip = sPolar.sel(x=slice(*[sxmin,sxmax]), y=slice(*[symin,symax]))

        sPolar.close()
        southpole.close()

        return self._interpUTM(southClip)

    def _getNorthPole(self, ds):

        # Polar Stereographic North (60N,0E) epsg:3995
        source_crs = 'epsg:3995' 
        target_crs = 'epsg:4326'
        latlon_to_npolar = pyproj.Transformer.from_crs(target_crs, source_crs)
        nblon, nblat = self._lon_lat_box([-180,180],[60,90], refinement=100)
        nbx, nby = latlon_to_npolar.transform(nblat, nblon)

        nxmin,nxmax = nbx.min(),nbx.max()
        nymin,nymax = nby.min(),nby.max()

        northpole = ds.sel(latitude=slice(*[35, 90]))
        northpole.rio.write_crs(target_crs, inplace=True)
        nPolar = northpole.rio.reproject(source_crs)
        nPolar = nPolar.sortby(nPolar.x)
        nPolar = nPolar.sortby(nPolar.y)
        northClip = nPolar.sel(x=slice(*[nxmin,nxmax]), y=slice(*[nymin,nymax]))

        nPolar.close()
        northpole.close()

        return self._interpUTM(northClip)

    def _cmptFlex(self, xrutm, regrid=1, shifted=False):

        if regrid>1:
            xcoord = xrutm.x.values.copy()
            ycoord = xrutm.y.values.copy()
            xrutm = xrutm.coarsen(y=regrid, x=regrid, 
                                    boundary='pad').mean()
        
        if shifted:
            cload = xrutm.shifted.values.copy()
        else:
            cload = xrutm.erodep.values.copy()

        z = xrutm.elevation.values.copy()
        shape = cload.shape
        dx = xrutm.x.values[1]-xrutm.x.values[0]
        dy = xrutm.y.values[1]-xrutm.y.values[0]

        flex = F2D()
        flex.Quiet = True

        flex.Method = "FD"  
        flex.PlateSolutionType = "vWC1994"  
        flex.Solver = "direct"  

        flex.g = 9.8  # acceleration due to gravity
        flex.E = 65e9  # Young's Modulus
        flex.nu = 0.25  # Poisson's Ratio
        flex.rho_m = 3300.0  # MantleDensity
        flex.rho_fill = 2300.0  # InfillMaterialDensity
        rho_water = 1000.0  # InfillMaterialDensity

        flex.Te = 50000.0 * np.ones(shape)  # Elastic thickness [m] -- scalar but may be an array

        flex.qs = cload * flex.rho_fill * flex.g
        r,c = np.where(z<0)
        flex.qs[r,c] = cload[r,c] * (flex.rho_fill-rho_water) * flex.g
        flex.dx = dx 
        flex.dy = dy

        # Boundary conditions can be:
        flex.BC_E = "Periodic"  # west boundary condition
        flex.BC_W = "Periodic"  # east boundary condition
        flex.BC_S = "Periodic"  # south boundary condition
        flex.BC_N = "Periodic"  # north boundary condition

        flex.initialize()
        flex.run()
        flex.finalize()
        xrutm['flex'] = (['y','x'],flex.w)

        if regrid>1:
            xrutm = xrutm.interp(x=xcoord, y=ycoord, method="nearest")
            xrutm.flex.rio.write_nodata(np.nan, inplace=True)
            xrutm.elevation.rio.write_nodata(np.nan, inplace=True)
            xrutm.erodep.rio.write_nodata(np.nan, inplace=True)
            xrutm.shiftz.rio.write_nodata(np.nan, inplace=True)
            xrutm.shifted.rio.write_nodata(np.nan, inplace=True)
            return xrutm.rio.interpolate_na() 

        return xrutm

    def _runFlexure(self, val_ds):

        # Hard-coded values could be changed in the future
        flexgrd = 2
        
        flexds = None
        tstart = process_time()
        tstep = process_time()
        val_ds = self._shiftValues(val_ds)

        # Split the global grids into tiles and project in Cartesian coordinates
        ##########################################################################
        lstbnds = self._makeBnds()
        listbnds = []
        if MPIsize == 1:
            listbnds = lstbnds.copy()
        elif MPIsize == 2:
            proclist = np.array_split(np.arange(len(lstbnds),dtype=int), MPIsize)
            listbnds = lstbnds[proclist[MPIrank][0]:proclist[MPIrank][-1]+1]
        else:
            proclist = np.array_split(np.arange(len(lstbnds),dtype=int), self.gflexproc-2)
            if MPIrank > 1 and MPIrank < self.gflexproc:
                listbnds = lstbnds[proclist[MPIrank-2][0]:proclist[MPIrank-2][-1]+1]

        list_clip_utm = []
        for k in range(len(listbnds)):
            dll = 10
            if listbnds[k][2] <= -30 or listbnds[k][2] >= 30:
                dll = 20
            if listbnds[k][2] <= -45 or listbnds[k][2] >= 45:
                dll = 30
            if listbnds[k][2] <= -60 or listbnds[k][2] >= 60:
                dll = 40
            tclip_utm = self._makeUTM(val_ds,listbnds[k],dll)
            list_clip_utm.append(tclip_utm)

        # Get polar grids
        if MPIrank == 0:
            north_utm = self._getNorthPole(val_ds)
        if MPIsize > 1:
            if MPIrank == 1:
                south_utm = self._getSouthPole(val_ds)
        else:
            south_utm = self._getSouthPole(val_ds)

        if MPIrank == 0 and self.verbose:
            print(
                "--- Build cartesian tiles (%0.02f seconds)"
                % (process_time() - tstep),
                flush=True,
            )

        # Compute flexural isostasy locally and reproject to lon/lat
        ##############################################################
        tstep = process_time()
        grd_flex = []
        bnd_flex = []
        grd_flex_edges = []
        bnd_flex_edges = []
        for k in range(len(listbnds)):
            if listbnds[k][0]>=-140 and listbnds[k][0]<140:
                bnd_flex.append(listbnds[k])
                grd_flex.append(self._cmptFlex(list_clip_utm[k],flexgrd,shifted=False))
            if listbnds[k][0]>=-40 and listbnds[k][0]<40:
                bnd_flex_edges.append(listbnds[k])
                grd_flex_edges.append(self._cmptFlex(list_clip_utm[k],flexgrd,shifted=True))

        flexGrids = []
        for k in range(len(grd_flex)):
            tmp = grd_flex[k]['flex']
            tmp = self._back2LonLat(tmp.to_dataset(),bnd_flex[k],self.res)
            flexGrids.append(tmp)  

        flexEdges = []
        for k in range(len(grd_flex_edges)):
            tmp = grd_flex_edges[k]['flex']
            tmp = self._back2LonLat(tmp.to_dataset(),bnd_flex_edges[k],self.res)
            flexEdges.append(tmp)      

        # Save local files on disk 
        if MPIrank > 0:
            if len(flexGrids) > 0:
                file_grids = self.outputDir+'/flexGrids_'+str(MPIrank)+'.pickle'
                with open(file_grids, 'wb') as file:
                    pickle.dump(flexGrids, file)
            if len(flexEdges) > 0:
                file_edges = self.outputDir+'/flexEdges_'+str(MPIrank)+'.pickle'
                with open(file_edges, 'wb') as file:
                    pickle.dump(flexEdges, file)

        # Run gflex on polar grids
        if MPIrank == 0:
            north_flex = self._cmptFlex(north_utm,flexgrd,shifted=False)
            tmp = self._makeLonLat(north_flex,[-180,180,75,90],self.res)
            flexNGrid = self._interpLonLat2(tmp,[-180,180,75,90],self.res)['flex'].to_dataset()
            flexNGrid.flex.rio.write_nodata(np.nan, inplace=True)
            flexNGrid = flexNGrid.rio.interpolate_na() 

        if MPIsize > 1:
            if MPIrank == 1:
                south_flex = self._cmptFlex(south_utm,flexgrd,shifted=False)
                tmp = self._makeLonLat(south_flex,[-180,180,-90,-75],self.res)
                flexSGrid = self._interpLonLat2(tmp,[-180,180,-90,-75],self.res)['flex'].to_dataset()
                flexSGrid.flex.rio.write_nodata(np.nan, inplace=True)
                flexSGrid = flexSGrid.rio.interpolate_na() 
                south_grid = self.outputDir+'/sthGrid.pickle'
                with open(south_grid, 'wb') as file:
                        pickle.dump(flexSGrid, file)
        else:
            south_flex = self._cmptFlex(south_utm,flexgrd,shifted=False)
            tmp = self._makeLonLat(south_flex,[-180,180,-90,-75],self.res)
            flexSGrid = self._interpLonLat2(tmp,[-180,180,-90,-75],self.res)['flex'].to_dataset()
            flexSGrid.flex.rio.write_nodata(np.nan, inplace=True)
            flexSGrid = flexSGrid.rio.interpolate_na() 

        MPIcomm.Barrier()
        if MPIrank == 0 and self.verbose:
                print(
                    "--- Perform flexural isostasy (%0.02f seconds)"
                    % (process_time() - tstep),
                    flush=True,
                )
        tstep = process_time()
        flexds = None
        if MPIrank == 0:

            # Combine the local grids into a single one
            flexds = None
            flexdse = None
            for k in range(1,len(flexGrids)):
                if k == 1:
                    flexds = flexGrids[0].combine_first(flexGrids[1])
                else:
                    flexds = flexds.combine_first(flexGrids[k])
            for k in range(1,len(flexEdges)):
                if k == 1:
                    flexdse = flexEdges[0].combine_first(flexEdges[1])
                else:
                    flexdse = flexdse.combine_first(flexEdges[k])  

            if MPIsize > 1:
                with open(self.outputDir+'/sthGrid.pickle', 'rb') as file:
                    flexSGrid = pickle.load(file)
            
            # Add the ones built on other processors
            for r in range(2,self.gflexproc):
                file_grids = self.outputDir+'/flexGrids_'+str(r)+'.pickle'
                if Path(file_grids).is_file():
                    rank_flexg = None
                    # Open the file in binary mode
                    with open(file_grids, 'rb') as file:
                        # Deserialize and retrieve the variable from the file
                        rank_flexg = pickle.load(file)
                        if flexds is None:
                            for k in range(1,len(rank_flexg)):
                                if k == 1:
                                    flexds = rank_flexg[0].combine_first(rank_flexg[1])
                                else:
                                    flexds = flexds.combine_first(rank_flexg[k])
                        else:
                            for k in range(0,len(rank_flexg)):
                                flexds = flexds.combine_first(rank_flexg[k])

                file_gridse = self.outputDir+'/flexEdges_'+str(r)+'.pickle'
                if Path(file_gridse).is_file():
                    rank_flexe = None
                    # Open the file in binary mode
                    with open(file_gridse, 'rb') as file:
                        # Deserialize and retrieve the variable from the file
                        rank_flexe = pickle.load(file)
                        if flexdse is None:
                            for k in range(1,len(rank_flexe)):
                                if k == 1:
                                    flexdse = rank_flexe[0].combine_first(rank_flexe[1])
                                else:
                                    flexdse = flexdse.combine_first(rank_flexe[k])
                        else:
                            for k in range(0,len(rank_flexe)):
                                flexdse = flexdse.combine_first(rank_flexe[k])
            
            flexds = flexds.combine_first(flexNGrid)
            flexds = flexds.combine_first(flexSGrid)
            flexds = self._getFlex(val_ds,flexds,flexdse) 

            for filename in glob.glob(self.outputDir+"/*.pickle"):
                os.remove(filename) 

            if self.verbose:
                print(
                    "Total isostatic flexural compensation runtime (%0.02f seconds)"
                    % (process_time() - tstart),
                    flush=True,
                )
        MPIcomm.Barrier()
        
        return flexds

    def globalFlexIso(self):

        if not self.gflexOn:
            return

        # Get elevations from the unstructured mesh structure 
        hl = self.hLocal.getArray().copy()
        newZ = np.zeros(self.mpoints, dtype=np.float64) - 1.0e8
        newZ[self.locIDs] = hl
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, newZ, op=MPI.MAX)
        
        # Get erosion deposition from the unstructured mesh structure
        edl = self.cumEDLocal.getArray().copy()
        prev_ed = self.cumEDFlex.getArray().copy()
        self.cumEDFlex.setArray(edl)
        gED = np.zeros(self.mpoints, dtype=np.float64) - 1.0e10
        gED[self.locIDs] = edl-prev_ed
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gED, op=MPI.MAX)

        # Interpolate to a lonlat regular grid
        if gED[self.gflex_ids].ndim == 2:
            erodep = (np.sum(self.gflex_wghts * gED[self.gflex_ids][:,:], axis=1)/self.gflex_sumwght) 
            elev = (np.sum(self.gflex_wghts * newZ[self.gflex_ids][:,:], axis=1)/self.gflex_sumwght)  
        else:
            erodep = (np.sum(self.gflex_wghts * gED[self.gflex_ids][:,:,0], axis=1)/self.gflex_sumwght)
            elev = (np.sum(self.gflex_wghts * newZ[self.gflex_ids][:,:,0], axis=1)/self.gflex_sumwght)  
        if len(self.gflex_oIDs) > 0:
            erodep[self.gflex_oIDs] = gED[self.gflex_ids[self.gflex_oIDs, 0]]
            elev[self.gflex_oIDs] = newZ[self.gflex_ids[self.gflex_oIDs, 0]]
        elev = np.reshape(elev, self.regshape)
        erodep = np.reshape(erodep, self.regshape)
        ds = xr.Dataset({
                'elevation': xr.DataArray(elev, coords=dict(latitude=self.lati, 
                                                longitude=self.long), 
                                dims=("latitude", "longitude")),
                'erodep': xr.DataArray(erodep, coords=dict(latitude=self.lati, 
                                                longitude=self.long), 
                                dims=("latitude", "longitude")),
                }
            )

        # Compute the flexural responses associated with the corresponding loads.
        flexds = self._runFlexure(ds)

        # Interpolate the calculate global thickness on the spherical mesh
        if MPIrank == 0:
            flexg = flexds.flex.values.flatten()
            uflex = np.sum(self.rflex_wghts * flexg[self.rflex_ids], axis=1) * self.rflex_denum
            if len(self.rflex_oIDs) > 0:
                uflex[self.rflex_oIDs] = flexg[self.rflex_ids[self.rflex_oIDs, 0]]
        else:
            uflex = None
        
        # Send flexural response globally
        globflex = MPI.COMM_WORLD.bcast(uflex, root=0)
        # Local flexural isostasy
        tmpFlex = globflex[self.locIDs]
        tmp = self.cumFlexL.getArray().copy()
        self.cumFlexL.setArray(tmp + tmpFlex)

        # Update elevation
        tmp = self.hLocal.getArray().copy()
        self.hLocal.setArray(tmp + tmpFlex)
        self.dm.localToGlobal(self.hLocal, self.hGlobal)

        # Store current erosion deposition values
        tmp = self.cumEDLocal.getArray().copy()
        self.cumEDFlex.setArray(tmp)

        self.nextFlex = self.tNow + self.gflexStep

        return