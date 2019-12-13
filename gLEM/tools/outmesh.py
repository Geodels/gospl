import gc
import h5py

import numpy as np
from mpi4py import MPI
import sys,petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
from time import clock

MPIrank = PETSc.COMM_WORLD.Get_rank()
MPIsize = PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD

try: range = xrange
except: pass

class WriteMesh(object):
    """
    Output model paramters using hdf5 library
    """
    def __init__(self):

        self.step = 0
        self.stratStep = 0
        self.file = 'gLEM'
        if MPIrank == 0 :
            self.create_OutputDir()

        # Sync the chosen output dir to all processors
        self.outputDir = MPIcomm.bcast(self.outputDir, root=0)

        if self.rStep > 0:
            self.step = self.rStep+1
            if self.strat > 0:
                n = int(self.tout/self.strat)
                self.stratStep = self.rStep*n
            self.tNow = self.tStart+self.rStep*self.tout
            self.rStart = self.tStart+self.rStep*self.tout
            self.saveTime = self.tNow+self.tout
            self.saveStrat = self.tNow+self.strat
        else:
            self.rStart = self.tStart

        return

    def create_OutputDir(self):
        """
        Create a directory to store outputs.
        """
        import os
        import glob
        import shutil

        # Get output directory
        if self.outputDir is not None:
            self.outputDir = os.getcwd()+'/'+self.outputDir
        else:
            self.outputDir = os.getcwd()+'/output'

        if self.rStep == 0:
            if self.makedir:
                if os.path.exists(self.outputDir):
                    self.outputDir += '_'+str(len(glob.glob(self.outputDir+str('*')))-1)
            else:
                if os.path.exists(self.outputDir):
                    shutil.rmtree(self.outputDir, ignore_errors=True)

            os.makedirs(self.outputDir)
            os.makedirs(self.outputDir+'/h5')
            os.makedirs(self.outputDir+'/xmf')

        return

    def outputStrat(self):
        """
        Saves mesh local stratigraphic information stored in the DMPlex to HDF5 file.
        """

        t = clock()
        h5file = self.outputDir+'/h5/stratal.'+str(self.stratStep)+'.p'+str(MPIrank)+'.h5'
        with h5py.File(h5file, "w") as f:
            f.create_dataset('elev',shape=(len(self.lcoords[:,0]),1), dtype='float32', compression='gzip')
            f["elev"][:,0] = self.hLocal.getArray()
            f.create_dataset('sea',shape=(1,1), dtype='float32', compression='gzip')
            f["sea"][:,0] = self.sealevel
            f.create_dataset('erodep',shape=(len(self.lcoords[:,0]),1), dtype='float32', compression='gzip')
            f["erodep"][:,0] = self.cumEDLocal.getArray()

        MPIcomm.Barrier()
        if MPIrank == 0 and self.verbose:
            print('Creating stratal outputfile (%0.02f seconds)'% (clock() - t))

        self.stratStep += 1

        return

    def outputMesh(self):
        """
        Saves mesh local information stored in the DMPlex to HDF5 file
        If the file already exists, it is overwritten.
        """

        t = clock()
        if self.step == 0:
            topology = self.outputDir+'/h5/topology.p'+str(MPIrank)+'.h5'
            with h5py.File(topology, "w") as f:
                f.create_dataset('coords',shape=(len(self.lcoords[:,0]),3), dtype='float32', compression='gzip')
                f["coords"][:,:] = self.lcoords
                f.create_dataset('cells',shape=(len(self.lcells[:,0]),3), dtype='int32', compression='gzip')
                f["cells"][:,:] = self.lcells+1
            self.elems = MPIcomm.gather(len(self.lcells[:,0]),root = 0)
            self.nodes = MPIcomm.gather(len(self.lcoords[:,0]),root = 0)

        h5file = self.outputDir+'/h5/'+self.file+'.'+str(self.step)+'.p'+str(MPIrank)+'.h5'
        with h5py.File(h5file, "w") as f:
            if self.stratStep == 0:
                f.create_dataset('elev',shape=(len(self.lcoords[:,0]),1), dtype='float32', compression='gzip')
                f["elev"][:,0] = self.hLocal.getArray()
                f.create_dataset('erodep',shape=(len(self.lcoords[:,0]),1), dtype='float32', compression='gzip')
                f["erodep"][:,0] = self.cumEDLocal.getArray()
            f.create_dataset('flowAcc',shape=(len(self.lcoords[:,0]),1), dtype='float32', compression='gzip')
            data = self.FAL.getArray().copy()
            data[data<=0.] = 1.
            f["flowAcc"][:,0] = data
            f.create_dataset('sedLoad',shape=(len(self.lcoords[:,0]),1), dtype='float32', compression='gzip')
            # data = self.vSedLocal.getArray().copy()
            # data[data<1.] = 1
            f["sedLoad"][:,0] = self.vSedLocal.getArray().copy()
            if self.uplift is not None:
                f.create_dataset('uplift',shape=(len(self.lcoords[:,0]),1), dtype='float32', compression='gzip')
                f["uplift"][:,0] = self.uplift
            if self.hdisp is not None:
                f.create_dataset('hdisp',shape=(len(self.lcoords[:,0]),3), dtype='float32', compression='gzip')
                f["hdisp"][:,:] = self.hdisp

            del data

        if MPIrank == 0:
            self._save_DMPlex_XMF()
            self._save_XDMF()

        MPIcomm.Barrier()
        if MPIrank == 0 and self.verbose:
            print('Creating outputfile (%0.02f seconds)'% (clock() - t))

        if MPIrank == 0:
            print('+++ Output Simulation Time: %0.02f years'% (self.tNow))

        self.step += 1
        gc.collect()

        return

    def readData(self):
        """
        For restarted simulations, this model reads the previous dataset.
        """
        import os

        if self.stratStep == 0:
            h5file = self.outputDir+'/h5/'+self.file+'.'+str(self.step)+'.p'+str(MPIrank)+'.h5'
        else:
            h5file = self.outputDir+'/h5/stratal.'+str(self.stratStep-1)+'.p'+str(MPIrank)+'.h5'

        if os.path.exists(h5file):
            hf = h5py.File(h5file, 'r')
        else:
            raise ValueError('Restart file is missing...')

        self.hLocal.setArray(np.array((hf['/elev'])))
        self.dm.localToGlobal(self.hLocal, self.hGlobal)
        self.cumEDLocal.setArray(np.array((hf['/erodep'])))
        self.dm.localToGlobal(self.cumEDLocal, self.cumED)

        self.elems = MPIcomm.gather(len(self.lcells[:,0]),root = 0)
        self.nodes = MPIcomm.gather(len(self.lcoords[:,0]),root = 0)

        hf.close()

        return

    def _save_DMPlex_XMF(self):
        """
        Saves mesh local information stored in the HDF5 to XMF file
        to visualise in Paraview.
        """

        xmf_file = self.outputDir+'/xmf/'+self.file+str(self.step)+'.xmf'

        f = open(xmf_file, 'w')

        # Header for xml file
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n')
        f.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
        f.write(' <Domain>\n')
        f.write('    <Grid GridType="Collection" CollectionType="Spatial">\n')
        f.write('      <Time Type="Single" Value="%0.02f"/>\n'%self.saveTime)

        for p in range(MPIsize):
            pfile = 'h5/'+str(self.file)+'.'+str(self.step)+'.p'+str(p)+'.h5'
            if self.stratStep == 0:
                afile  = 'h5/'+str(self.file)+'.'+str(self.step)+'.p'+str(p)+'.h5'
            else:
                afile = 'h5/stratal.'+str(self.stratStep-1)+'.p'+str(p)+'.h5'
            tfile = 'h5/topology.p'+str(p)+'.h5'
            f.write('      <Grid Name="Block.%s">\n' %(str(p)))
            f.write('         <Topology Type="Triangle" NumberOfElements="%d" BaseOffset="1">\n'%self.elems[p])
            f.write('          <DataItem Format="HDF" DataType="Int" ')
            f.write('Dimensions="%d 3">%s:/cells</DataItem>\n'%(self.elems[p],tfile))
            f.write('         </Topology>\n')

            f.write('         <Geometry Type="XYZ">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 3">%s:/coords</DataItem>\n'%(self.nodes[p],tfile))
            f.write('         </Geometry>\n')

            f.write('         <Attribute Type="Scalar" Center="Node" Name="Z">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/elev</DataItem>\n'%(self.nodes[p],afile))
            f.write('         </Attribute>\n')

            f.write('         <Attribute Type="Scalar" Center="Node" Name="ED">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/erodep</DataItem>\n'%(self.nodes[p],afile))
            f.write('         </Attribute>\n')

            f.write('         <Attribute Type="Scalar" Center="Node" Name="FA">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/flowAcc</DataItem>\n'%(self.nodes[p],pfile))
            f.write('         </Attribute>\n')

            f.write('         <Attribute Type="Scalar" Center="Node" Name="SL">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/sedLoad</DataItem>\n'%(self.nodes[p],pfile))
            f.write('         </Attribute>\n')
            if self.hdisp is not None:
                f.write('         <Attribute Type="Vector" Center="Node" Name="hTec">\n')
                f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
                f.write('Dimensions="%d 3">%s:/hdisp</DataItem>\n'%(self.nodes[p],pfile))
                f.write('         </Attribute>\n')
            if self.uplift is not None:
                f.write('         <Attribute Type="Scalar" Center="Node" Name="vTec">\n')
                f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
                f.write('Dimensions="%d 1">%s:/uplift</DataItem>\n'%(self.nodes[p],pfile))
                f.write('         </Attribute>\n')

            f.write('         <Attribute Type="Scalar" Center="Node" Name="sea">\n')
            f.write('          <DataItem ItemType="Function" Function="$0 * 0.00000000001 + %f" Dimensions="%d 1">\n'%(self.sealevel,self.nodes[p]))
            f.write('           <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/erodep</DataItem>\n'%(self.nodes[p],pfile))
            f.write('          </DataItem>\n')
            f.write('         </Attribute>\n')

            f.write('      </Grid>\n')

        f.write('    </Grid>\n')
        f.write(' </Domain>\n')
        f.write('</Xdmf>\n')
        f.close()

        return

    def _save_XDMF(self):
        """
        This function writes the XDmF file which is calling the XmF file.
        """

        xdmf_file = self.outputDir+'/'+self.file+'.xdmf'
        f= open(xdmf_file,'w')

        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n')
        f.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
        f.write(' <Domain>\n')
        f.write('    <Grid GridType="Collection" CollectionType="Temporal">\n')

        for s in range(self.step+1):
            xmf_file = 'xmf/'+self.file+str(s)+'.xmf'
            f.write('      <xi:include href="%s" xpointer="xpointer(//Xdmf/Domain/Grid)"/>\n' %xmf_file)

        f.write('    </Grid>\n')
        f.write(' </Domain>\n')
        f.write('</Xdmf>\n')
        f.close()

        return

    # def outputMesh(self):
    #     """
    #     Saves mesh local information stored in the DMPlex to HDF5 file
    #     If the file already exists, it is overwritten.
    #     """
    #
    #     t = clock()
    #     if self.step == 0:
    #         topology = self.outputDir+'/h5/topology.p'+str(MPIrank)+'.h5'
    #         with h5py.File(topology, "w") as f:
    #             f.create_dataset('coords',shape=(len(self.lcoords[:,0]),3), dtype='float32', compression='gzip')
    #             f["coords"][:,:] = self.lcoords
    #             f.create_dataset('cells',shape=(len(self.lcells[:,0]),3), dtype='int32', compression='gzip')
    #             f["cells"][:,:] = self.lcells+1
    #         self.elems = MPIcomm.gather(len(self.lcells[:,0]),root = 0)
    #         self.nodes = MPIcomm.gather(len(self.lcoords[:,0]),root = 0)
    #
    #     h5file = self.outputDir+'/h5/'+self.file+'.'+str(self.step)+'.p'+str(MPIrank)+'.h5'
    #     with h5py.File(h5file, "w") as f:
    #         f.create_dataset('elev',shape=(len(self.lcoords[:,0]),1), dtype='float32', compression='gzip')
    #         f["elev"][:,0] = self.hLocal.getArray()
    #         f.create_dataset('flowAcc',shape=(len(self.lcoords[:,0]),1), dtype='float32', compression='gzip')
    #         data = self.FAL.getArray().copy()
    #         data[data<=0.] = 1.
    #         f["flowAcc"][:,0] = data
    #         f.create_dataset('erodep',shape=(len(self.lcoords[:,0]),1), dtype='float32', compression='gzip')
    #         f["erodep"][:,0] = self.cumEDLocal.getArray()
    #         f.create_dataset('sedLoad',shape=(len(self.lcoords[:,0]),1), dtype='float32', compression='gzip')
    #         # data = self.vSedLocal.getArray().copy()
    #         # data[data<1.] = 1
    #         f["sedLoad"][:,0] = self.vSedLocal.getArray().copy()
    #         del data
    #
    #     if MPIrank == 0:
    #         self._save_DMPlex_XMF()
    #         self._save_XDMF()
    #
    #     MPIcomm.Barrier()
    #     if MPIrank == 0 and self.verbose:
    #         print('Creating outputfile (%0.02f seconds)'% (clock() - t))
    #
    #     if MPIrank == 0:
    #         print('+++ Output Simulation Time: %0.02f years'% (self.tNow))
    #
    #     self.step += 1
    #     gc.collect()
    #
    #     return
