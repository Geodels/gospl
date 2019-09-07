import numpy as np
import pandas as pd
import ruamel.yaml as yaml
from scipy.interpolate import interp1d

from mpi4py import MPI
import sys,petsc4py
petsc4py.init(sys.argv)

from petsc4py import PETSc

MPIrank = PETSc.COMM_WORLD.Get_rank()
MPIsize = PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD

class ReadYaml(object):
    """
    Parsing YAML input file

    Args:
        filename: input filename (.yml YAML file)
    """
    def __init__(self, filename):

        # Check input file exists
        try:
            with open(filename) as finput:
                pass
        except IOError as exc:
            print("Unable to open file: ",filename)
            raise IOError('The input file is not found...')

        # Open YAML file
        with open(filename, 'r') as finput:
            self.input = yaml.load(finput, Loader=yaml.Loader)

        if MPIrank == 0 and 'name' in self.input.keys() and self.verbose:
            print("The following model will be run:     {}".format(self.input['name']))

        # Read simulation parameters
        self._readDomain()
        self._readTime()
        self._readSPL()
        self._readHillslope()
        self._readSealevel()
        self._readOut()

        self.tNow = self.tStart
        self.saveTime = self.tNow

        return

    def _readDomain(self):
        """
        Read domain definition, boundary conditions and flow direction parameters.
        """

        try:
            domainDict = self.input['domain']
        except KeyError as exc:
            print("Key 'domain' is required and is missing in the input file!")
            raise KeyError('Key domain is required in the input file!')

        try:
            self.radius = domainDict['radius']
        except KeyError as exc:
            self.radius = 6378137.0

        try:
            self.flowDir = domainDict['flowdir']
        except KeyError as exc:
            self.flowDir = 6

        try:
            self.reflevel = domainDict['refinement']
        except KeyError as exc:
            self.reflevel = 9

        try:
            meshFile = domainDict['npdata']
        except KeyError as exc:
            print("Key 'npdata' is required and is missing in the 'domain' declaration!")
            raise KeyError('Compressed numpy dataset definition is not defined!')

        self.meshFile = meshFile+str(self.reflevel)+'.npz'
        try:
            with open(self.meshFile) as meshfile:
                pass
        except IOError as exc:
            print("Unable to open numpy dataset: ",self.meshFile)
            raise IOError('The numpy dataset is not found...')

        return

    def _readTime(self):
        """
        Read simulation time declaration.
        """

        try:
            timeDict = self.input['time']
        except KeyError as exc:
            print("Key 'time' is required and is missing in the input file!")
            raise KeyError('Key time is required in the input file!')

        try:
            self.tStart = timeDict['start']
        except KeyError as exc:
            print("Key 'start' is required and is missing in the 'time' declaration!")
            raise KeyError('Simulation start time needs to be declared.')

        try:
            self.tEnd = timeDict['end']
        except KeyError as exc:
            print("Key 'end' is required and is missing in the 'time' declaration!")
            raise KeyError('Simulation end time needs to be declared.')

        try:
            self.dt = timeDict['dt']
        except KeyError as exc:
            print("Key 'dt' is required and is missing in the 'time' declaration!")
            raise KeyError('Simulation discretisation time step needs to be declared.')

        try:
            self.tout = timeDict['tout']
        except KeyError as exc:
            self.tout = self.tEnd-self.tStart
            print("Output time interval 'tout' has been set to {} years".format(self.tout))

        if self.tEnd <= self.tStart:
            raise ValueError('Simulation end/start times do not make any sense!')

        if self.tout < self.dt:
            self.tout = self.dt
            print("Output time interval was changed to {} years to match the time step dt".format(self.dt))

        return

    def _readSPL(self):
        """
        Read surface processes bedrock parameters.
        """

        try:
            splDict = self.input['spl']
            try:
                self.K = splDict['K']
            except KeyError as exc:
                print("When using the Surface Process Model definition of coefficient Kb is required.")
                raise ValueError('Surface Process Model: Kb coefficient not found.')
            try:
                self.frac_fine = splDict['Ff']
            except KeyError as exc:
                self.frac_fine = 0.0
            try:
                # `wgth` is the percentage of upstream sediment flux that will be deposited on each cell...
                self.wgth = splDict['wgth']
                if self.wgth>=1.0:
                    self.wgth = 0.999
            except KeyError as exc:
                self.wgth = 0.

        except KeyError as exc:
            self.K = 1.e-12
            self.wgth = 0.
            self.frac_fine = 0.0

        return

    def _readHillslope(self):
        """
        Read hillslope parameters.
        # """

        try:
            hillDict = self.input['diffusion']
            try:
                self.Cd = hillDict['hillslopeK']
            except KeyError as exc:
                print("When declaring diffusion processes, the coefficient hillslopeK is required.")
                raise ValueError('Hillslope: Cd coefficient not found.')
            try:
                self.sedimentK = hillDict['sedimentK']
            except KeyError as exc:
                self.sedimentK = 10.
        except KeyError as exc:
            self.Cd = 0.
            self.sedimentK = 10.

        return

    def _readSealevel(self):
        """
        Define sealevel evolution.
        """

        seafile = None
        seafunction = None
        sealevel = 0.
        self.seafunction = None
        try:
            seaDict = self.input['sea']
            try:
                sealevel = seaDict['position']
                try:
                    seafile = seaDict['curve']
                except KeyError as exc:
                    seafile = None
            except KeyError as exc:
                try:
                    seafile = seaDict['curve']
                except KeyError as exc:
                    seafile = None
        except KeyError as exc:
            sealevel = 0.

        if seafile is not None:
            try:
                with open(seafile) as fsea:
                    try:
                        seadata = pd.read_csv(seafile, sep=r',', engine='c',
                                                  header=None, na_filter=False,
                                                  dtype=np.float, low_memory=False)
                        pass
                    except ValueError as exc:
                        try:
                            seadata = pd.read_csv(seafile, sep=r'\s+',
                                                      engine='c', header=None,
                                                      na_filter=False, dtype=np.float,
                                                      low_memory=False)
                            pass
                        except ValueError as exc:
                            print("The sea-level file is not well formed: it should be comma or tab separated")
                            raise ValueError('Wrong formating of sea-level file.')
            except IOError as exc:
                print("Unable to open file: ",seafile)
                raise IOError('The sealevel file is not found...')

            seadata[1] += sealevel
            if seadata[0].min() > self.tStart:
                tmpS = []
                tmpS.insert(0, {0: self.tStart, 1: seadata[1].iloc[0]})
                seadata = pd.concat([pd.DataFrame(tmpS), seadata], ignore_index=True)
            if seadata[0].max() < self.tEnd:
                tmpE = []
                tmpE.insert(0, {0: self.tEnd, 1: seadata[1].iloc[-1]})
                seadata = pd.concat([seadata,pd.DataFrame(tmpE)], ignore_index=True)
            self.seafunction = interp1d(seadata[0], seadata[1]+sealevel, kind='linear')
        else:
            year = np.linspace(self.tStart, self.tEnd+self.dt, num=11, endpoint=True)
            seaval = np.full(len(year),sealevel)
            self.seafunction = interp1d(year, seaval, kind='linear')

        return

    def _readOut(self):
        """
        Parse output directory.
        """

        try:
            outDict = self.input['output']
            try:
                self.outputDir = outDict['dir']
            except KeyError as exc:
                self.outputDir = 'output'
            try:
                self.makedir = outDict['makedir']
            except KeyError as exc:
                self.makedir = True
        except KeyError as exc:
            self.outputDir = 'output'
            self.makedir = True

        return
