import gc
import glob
import h5py
import math
import fileinput
import numpy as np
import pandas as pd
from mpi4py import MPI
from scipy import sparse
from scipy import spatial

from scipy import ndimage
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
import scipy.optimize
from scipy.linalg import norm
from scipy import sum, average
from scipy.stats import pearsonr
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

import sys,petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
from time import clock

from gLEM._fortran import fillPIT
from gLEM._fortran import MFDreceivers
from gLEM._fortran import setHillslopeCoeff
from gLEM._fortran import setDiffusionCoeff

MPIrank = PETSc.COMM_WORLD.Get_rank()
MPIsize = PETSc.COMM_WORLD.Get_size()
MPIcomm = PETSc.COMM_WORLD

try: range = xrange
except: pass

class PFit(object):
    """
    Constraining model outputs based on known final paleotopography map.
    """
    def __init__(self):

        self.runNb = 0
        self.simtime = self.tEnd - self.tStart
        self.accuracy = []
        self.similarity = []
        self.improvement = []

        # Get paleomap
        loadData = np.load(self.paleoTopo)
        gZ = loadData['z']
        if abs(len(gZ)-self.gpoints) > 0:
            print('Forcing paleotopography file size differs from loaded one: %s',self.paleoTopo)
            raise RuntimeError("Error in forcing paleotopography file!")

        self.paleoGlobal = self.dm.createGlobalVector()
        self.paleoLocal = self.dm.createLocalVector()
        self.diffZ = self.paleoLocal.duplicate()

        # From global values to local ones...
        self.paleoLocal.setArray(gZ[self.glIDs])
        self.dm.localToGlobal(self.paleoLocal, self.paleoGlobal)

        # Get longitude/latitude in degrees
        loadData = np.load(self.lonlat)
        icoords = loadData['v']
        icoords[0,:] = icoords[0,:]*360./3601.-180.
        icoords[1,:] = icoords[1,:]*180./1801.-90.

        # Target grid to interpolate to
        ilon = np.arange(-180.,180.,1.)
        ilat = np.arange(-90.,90.,1.)
        self.ilon,self.ilat = np.meshgrid(ilon,ilat)

        # KDtree interpolation parameters
        LL = np.zeros((self.ilon.size,2))
        LL[:,0] = self.ilon.ravel()
        LL[:,1] = self.ilat.ravel()
        self.fixtree = spatial.cKDTree(self.gCoords,leafsize=10)
        tree = spatial.cKDTree(icoords.T,leafsize=10)
        dist, self.ll_id = tree.query(LL, k=3)

        self.ll_wght = np.divide(np.ones(dist.shape), dist**2, out=np.zeros_like(dist), where=dist!=0)
        self.ll_oid = np.where(dist[:,0] == 0)[0]

        # Map paleotopography onto the regular grid
        self.paleoZ = np.sum(self.ll_wght*gZ[self.ll_id],axis=1)/np.sum(self.ll_wght, axis=1)
        if len(self.ll_oid)>0:
            self.paleoZ[self.ll_oid] = gZ[self.ll_id[self.ll_oid,0]]
        self.paleoZ = self.paleoZ.reshape(self.ilon.shape)

        self.paleoLand = np.where(self.paleoZ>=0.)
        self.paleoOcean = np.where(self.paleoZ<=-1000.)
        self.paleoShelf = np.where(np.logical_and(self.paleoZ<0.,self.paleoZ>-1000.))

        nbdisp = len(self.tecdata)
        self.init_zdisp = []
        for k in range(nbdisp):
            self.init_zdisp.append(self.tecdata.iloc[k,2])

        del loadData, gZ, ilon, ilat, tree, LL, dist, icoords
        gc.collect()

        return

    def comparePaleomaps(self,verbose=False,disk=False):
        """
        Compare local differences between modelled and reference paleomaps
        """

        # Compute the difference between paleo and simulated elevation:
        #   + positive: paleo elevation is above simulated one
        #   + negative: paleo elevation is below simulated one
        self.diffZ.waxpy(-1.0,self.hLocal,self.paleoLocal)

        # Define differences to reference paleomap in vertical displacement rates (m/yr)
        tmp = self.diffZ.getArray().copy()/self.simtime
        diff = np.zeros(self.gpoints)
        diff = tmp[self.lgIDs]
        diff[self.outIDs] = -1.e8
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, diff, op=MPI.MAX)

        # Map spherical mesh differences on the 1 degree lon/lat regular grid
        self.ll_diff = np.sum(self.ll_wght*diff[self.ll_id],axis=1)/np.sum(self.ll_wght, axis=1)
        if len(self.ll_oid)>0:
            self.ll_diff[self.ll_oid] = diff[self.ll_id[self.ll_oid,0]]
        self.ll_diff = self.ll_diff.reshape(self.ilon.shape)

        del tmp
        gc.collect()

        self._accuracyScores()

        # Improvement compared to previous model
        self.force_diff = None
        if self.runNb == 0:
            self.updated_diff = diff.copy()
        else:
            self._improvementModels()
            if self.accuracy[-2].values[0][4] >= self.cvglimit and self.runNb <= self.paleostep:
            # if int(self.improvement[-1].values[0][1]) >= self.cvglimit and self.runNb <= self.paleostep:
                self.updated_diff += diff
            else:
                self.force_diff = diff.copy()

        if verbose and MPIrank==0:
            print('Model run: ',self.runNb)
            print(self.similarity[-1])
            print('')
            print(self.accuracy[-1])
            print('')
            if self.runNb > 0:
                print(self.improvement[-1])
                print('')

        if disk:
            simDF = pd.concat(self.similarity)
            simDF.to_pickle('similarity')
            accDF = pd.concat(self.accuracy)
            accDF.to_pickle('accuracy')
            impDF = pd.concat(self.improvement)
            impDF.to_pickle('improvement')
            # if len(self.similarity)>2:
            #     simDF = pd.concat(self.similarity[:-2])
            #     simDF.to_pickle('similarity')
            #     accDF = pd.concat(self.accuracy[:-2])
            #     accDF.to_pickle('accuracy')
            #     impDF = pd.concat(self.improvement[:-2])
            #     impDF.to_pickle('improvement')
            # else:
            #     simDF = pd.concat(self.similarity)
            #     simDF.to_pickle('similarity')
            #     accDF = pd.concat(self.accuracy)
            #     accDF.to_pickle('accuracy')
            #     impDF = pd.concat(self.improvement)
            #     impDF.to_pickle('improvement')

        return

    def forcePaleomap(self):
        """
        Force the paleomap into the simulated final elevation
        """

        self.diffZ.waxpy(-1.0,self.hLocal,self.paleoLocal)
        self.paleoLocal.copy(result=self.hLocal)
        self.dm.localToGlobal(self.hLocal, self.hGlobal)

        return

    def backwardAdvection(self, verbose=False):
        """
        Advect vertical displacement missmatch backward in time
        """

        outfile = []
        nbdisp = len(self.tecdata)
        for k in range(nbdisp):
            outf = self.paleoDir+'/vdisp_r'+str(self.runNb)+'_s'+str(k)+'.npz'
            outfile.append(outf)

        if MPIrank == 0:
            advectXYZ = self.gCoords.copy()
            p = nbdisp-1
            for k in range(nbdisp):
                # Advect backward in time to previous coordinates
                if p < nbdisp-1:
                    dt = self.tecdata.iloc[p+1,0]-self.tecdata.iloc[p,0]
                else:
                    dt = self.tEnd-self.tecdata.iloc[p,0]

                if k==0 :
                    # Reading the last advection velocity
                    vel3d = np.load(self.tecdata.iloc[p,1])
                    advectXYZ -= vel3d['xyz']*dt
                else :
                    advectXYZ -= nadvect*dt

                # Interpolate from backward advected coordinates
                mvtree = spatial.cKDTree(advectXYZ,leafsize=10)
                # distances, indices = mvtree.query(self.gCoords, k=1)
                # nzdisp = self.updated_diff[indices]
                distances, indices = mvtree.query(self.gCoords, k=3)
                # Inverse weighting distance...
                weights = np.divide(np.ones(distances.shape), distances**2, out=np.zeros_like(distances), where=distances!=0)
                onIDs = np.where(distances[:,0] == 0)[0]
                nzdisp = np.sum(weights*self.updated_diff[indices],axis=1)/np.sum(weights, axis=1)
                if len(onIDs)>0:
                    nzdisp[onIDs] = self.updated_diff[indices[onIDs,0]]

                if self.init_zdisp[p] != 'empty':
                    # Previous vertical forcing file
                    disp0 = np.load(self.init_zdisp[p])
                    nzdisp += disp0['z']

                # Write vertical displacement at previous time:
                if p == nbdisp-1 and self.force_diff is not None:
                    np.savez_compressed(outfile[k],z=nzdisp+self.force_diff*self.simtime/dt)
                else:
                    np.savez_compressed(outfile[k], z=nzdisp)
                if verbose:
                    print("Processing {}".format(outfile[k]+".npz"))

                # Reading the previous advection velocity
                if p == 0:
                    break

                vel3d = np.load(self.tecdata.iloc[p-1,1])
                ndisps = vel3d['xyz']

                distances, indices = self.fixtree.query(advectXYZ, k=3)
                nadvect = np.zeros(advectXYZ.shape)

                # Inverse weighting distance...
                weights = np.divide(np.ones(distances.shape), distances**2, out=np.zeros_like(distances), where=distances!=0)

                onIDs = np.where(distances[:,0] == 0)[0]
                nadvect[:,0] = np.sum(weights*ndisps[indices,0],axis=1)/np.sum(weights, axis=1)
                nadvect[:,1] = np.sum(weights*ndisps[indices,1],axis=1)/np.sum(weights, axis=1)
                nadvect[:,2] = np.sum(weights*ndisps[indices,2],axis=1)/np.sum(weights, axis=1)

                if len(onIDs)>0:
                    nadvect[onIDs,0] = ndisps[indices[onIDs,0],0]
                    nadvect[onIDs,1] = ndisps[indices[onIDs,0],1]
                    nadvect[onIDs,2] = ndisps[indices[onIDs,0],2]

                p -= 1

        MPIcomm.Barrier()

        # Update vertical displacement files
        p = 0
        for k in reversed(range(nbdisp)):
            self.tecdata.iloc[p,2] = outfile[k]
            p += 1

        return

    def _computeNorms(self, map1, map2):
        """
        Manhattan distance is a distance metric between two points in a N dimensional vector space.
        It is the sum of the lengths of the projections of the line segment between the points onto the
        coordinate axes. In simple terms, it is the sum of absolute difference between the measures in
        all dimensions of two points.

        :arg map1: observation map
        :arg map2: simulation map
        """

        # Normalize
        rng = map1.max()-map1.min()
        amin = map1.min()
        map1 = (map1-amin)*255/rng
        rng = map2.max()-map2.min()
        amin = map2.min()
        map2 = (map2-amin)*255/rng

        return sum(abs(map2 - map1))

    def _similarity(self, X, Y=None, *, normalise=True, demean=True):
        """
        Compute similarity between the columns of one or two matrices.
         - Covariance: normalise=False, demean=True - https://en.wikipedia.org/wiki/Covariance
         - Corrcoef: normalise=True, demean=True - https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient
         - Dot product: normalise=False, demean=False
         - Cosine similarity: normalise=True, demean=False - https://en.wikipedia.org/wiki/Cosine_similarity N.B. also known as the congruence coefficient https://en.wikipedia.org/wiki/Congruence_coefficient
        """

        eps = 1.e-5
        if Y is None:
            if X.ndim != 2:
                raise ValueError("X must be 2D!")
            Y = X

        if X.ndim != 2 or Y.ndim != 2 or X.shape[0] != Y.shape[0]:
            raise ValueError("X and Y must be 2D with the same first dimension!")

        if demean:
            X = X - np.mean(X, axis=0)
            Y = Y - np.mean(Y, axis=0)

        if normalise:
            # Set variances to unity
            x = np.sqrt(np.sum(X**2, axis=0)); x[x < eps] = 1.0
            y = np.sqrt(np.sum(Y**2, axis=0)); y[y < eps] = 1.0
            X = X / x; Y = Y / y
        else:
            # Or just divide by no. of observations to make an expectation
            X = X / math.sqrt(X.shape[0]); Y = Y / math.sqrt(Y.shape[0])

        return X.T @ Y

    def _rvCoefficient(self, X, Y, *, normalise=True, demean=True):
        """
        RV coefficient (and related variants) between 2D matrices, columnwise
        Check: https://en.wikipedia.org/wiki/RV_coefficient

        RV normally defined in terms of corrcoefs, but any of the above similarity metrics will work
        """

        # Calculate correlations
        # N.B. trace(Xt @ Y) = sum(X * Y)
        Sxy = self._similarity(X, Y, normalise=normalise, demean=demean)
        c_xy = np.sum(Sxy ** 2)
        Sxx = self._similarity(X, X, normalise=normalise, demean=demean)
        c_xx = np.sum(Sxx ** 2)
        Syy = self._similarity(Y, Y, normalise=normalise, demean=demean)
        c_yy = np.sum(Syy ** 2)

        # And put together
        rv =  c_xy / math.sqrt(c_xx * c_yy)

        return rv

    def _concordance_correlation_coefficient(self, y_true, y_pred, sample_weight=None, multioutput='uniform_average'):
        """
        Concordance correlation coefficient.

        The concordance correlation coefficient is a measure of inter-rater agreement.
        It measures the deviation of the relationship between predicted and true values from the 45 degree angle.

        Read more: https://en.wikipedia.org/wiki/Concordance_correlation_coefficient
        Original paper: Lawrence, I., and Kuei Lin. "A concordance correlation coefficient to evaluate reproducibility." Biometrics (1989): 255-268.

        - y_true : array-like of shape = (n_samples) or (n_samples, n_outputs). Ground truth (correct) target values.
        - y_pred : array-like of shape = (n_samples) or (n_samples, n_outputs). Estimated target values.

        Returns
        - loss : A float in the range [-1,1]. A value of 1 indicates perfect agreement between the true and the predicted values.
        """

        cor = np.corrcoef(y_true,y_pred)[0][1]

        mean_true = np.mean(y_true)
        mean_pred = np.mean(y_pred)

        var_true = np.var(y_true)
        var_pred = np.var(y_pred)

        sd_true = np.std(y_true)
        sd_pred = np.std(y_pred)

        numerator = 2*cor*sd_true*sd_pred

        denominator = var_true+var_pred+(mean_true-mean_pred)**2

        return numerator/denominator

    def _nonMatching(self, diff):
        """
        Find non matching nodes.
        """

        return len(np.where(abs(diff)>self.erange)[0])/diff.size

    def _improvementModels(self):
        """
        Check improvement between consecutive simulations based on RMSE score
        """

        rmse1 = self.accuracy[-2]['RMSE']
        rmse2 = self.accuracy[-1]['RMSE']
        regions = ['global','lands','shelf','ocean']
        RI = (rmse1-rmse2)*100./rmse1
        fRI = pd.DataFrame({'region':regions, 'relative improvement': RI})
        self.improvement.append(fRI)

        return

    def _accuracyScores(self):

        obs = self.paleoZ.copy()
        sim = obs - self.ll_diff*self.simtime

        regions = ['global','lands','shelf','ocean']

        L1 = []
        CV = []
        CR = []
        DP = []
        CS = []
        P = []
        R2 = []
        MAE = []
        RMSE = []
        LCC = []
        NN = []

        # Global
        CV.append(self._rvCoefficient(obs, sim, normalise=False, demean=True))
        CR.append(self._rvCoefficient(obs, sim, normalise=True, demean=True))
        DP.append(self._rvCoefficient(obs, sim, normalise=False, demean=False))
        CS.append(self._rvCoefficient(obs, sim, normalise=True, demean=False))

        L1.append(self._computeNorms(obs,sim))
        P.append(pearsonr(obs.ravel(), sim.ravel())[0])
        R2.append(r2_score(obs.ravel(), sim.ravel(),multioutput='variance_weighted'))
        MAE.append(mean_absolute_error(obs.ravel(), sim.ravel(), multioutput='uniform_average'))
        RMSE.append(np.sqrt(mean_squared_error(obs.ravel(), sim.ravel())))
        LCC.append(self._concordance_correlation_coefficient(obs.ravel(), sim.ravel()))
        NN.append(self._nonMatching(sim-obs))

        # Land
        Lobs = obs.copy()
        Lsim = sim.copy()
        Lobs[self.paleoZ<0] = 0.
        Lsim[self.paleoZ<0] = 0.
        CV.append(self._rvCoefficient(Lobs, Lsim, normalise=False, demean=True))
        CR.append(self._rvCoefficient(Lobs, Lsim, normalise=True, demean=True))
        DP.append(self._rvCoefficient(Lobs, Lsim, normalise=False, demean=False))
        CS.append(self._rvCoefficient(Lobs, Lsim, normalise=True, demean=False))

        land_obs = obs[self.paleoLand].ravel()
        land_sim = sim[self.paleoLand].ravel()
        L1.append(self._computeNorms(land_obs,land_sim))
        P.append(pearsonr(land_obs.ravel(), land_sim.ravel())[0])
        R2.append(r2_score(land_obs, land_sim,multioutput='variance_weighted'))
        MAE.append(mean_absolute_error(land_obs, land_sim, multioutput='uniform_average'))
        RMSE.append(np.sqrt(mean_squared_error(land_obs, land_sim)))
        LCC.append(self._concordance_correlation_coefficient(land_obs, land_sim))
        NN.append(self._nonMatching(land_sim-land_obs))

        # Shelf
        Sobs = obs.copy()
        Ssim = sim.copy()
        Sobs[self.paleoZ>0] = 0.
        Ssim[self.paleoZ>0] = 0.
        Sobs[self.paleoZ<-1000.] = 0.
        Ssim[self.paleoZ<-1000.] = 0.
        CV.append(self._rvCoefficient(Sobs, Ssim, normalise=False, demean=True))
        CR.append(self._rvCoefficient(Sobs, Ssim, normalise=True, demean=True))
        DP.append(self._rvCoefficient(Sobs, Ssim, normalise=False, demean=False))
        CS.append(self._rvCoefficient(Sobs, Ssim, normalise=True, demean=False))

        shelf_obs = obs[self.paleoShelf].ravel()
        shelf_sim = sim[self.paleoShelf].ravel()
        L1.append(self._computeNorms(shelf_obs,shelf_sim))
        P.append(pearsonr(shelf_obs, shelf_sim)[0])
        R2.append(r2_score(shelf_obs, shelf_sim,multioutput='variance_weighted'))
        MAE.append(mean_absolute_error(shelf_obs, shelf_sim, multioutput='uniform_average'))
        RMSE.append(np.sqrt(mean_squared_error(shelf_obs, shelf_sim)))
        LCC.append(self._concordance_correlation_coefficient(shelf_obs, shelf_sim))
        NN.append(self._nonMatching(shelf_sim-shelf_obs))

        # Ocean
        Oobs = obs.copy()
        Osim = sim.copy()
        Lobs[self.paleoZ>=-1000.] = 0.
        Lsim[self.paleoZ>=-1000.] = 0.
        CV.append(self._rvCoefficient(Oobs, Osim, normalise=False, demean=True))
        CR.append(self._rvCoefficient(Oobs, Osim, normalise=True, demean=True))
        DP.append(self._rvCoefficient(Oobs, Osim, normalise=False, demean=False))
        CS.append(self._rvCoefficient(Oobs, Osim, normalise=True, demean=False))

        ocean_obs = obs[self.paleoOcean].ravel()
        ocean_sim = sim[self.paleoOcean].ravel()
        L1.append(self._computeNorms(ocean_obs,ocean_sim))
        P.append(pearsonr(ocean_obs.ravel(), ocean_sim.ravel())[0])
        R2.append(r2_score(ocean_obs.ravel(), ocean_sim.ravel(),multioutput='variance_weighted'))
        MAE.append(mean_absolute_error(ocean_obs.ravel(), ocean_sim.ravel(), multioutput='uniform_average'))
        RMSE.append(np.sqrt(mean_squared_error(ocean_obs.ravel(), ocean_sim.ravel())))
        LCC.append(self._concordance_correlation_coefficient(ocean_obs.ravel(), ocean_sim.ravel()))
        NN.append(self._nonMatching(ocean_sim-ocean_obs))

        self.similarity.append(pd.DataFrame({'region':regions, 'covariance': CV, 'corrcoef': CR,
                                 'dotproduct': DP, 'cosine': CS,
                                 'pearson': P, 'L1': L1}))

        self.accuracy.append(pd.DataFrame({'region':regions,  'RMSE': RMSE, 'R2': R2,
                         'MAE': MAE,  'LCC': LCC, 'nonmatching': NN}))

        return
