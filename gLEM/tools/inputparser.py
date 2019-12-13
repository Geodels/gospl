import os, errno
import numpy as np
import pandas as pd
import ruamel.yaml as yaml
from operator import itemgetter
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
        self._readTectonic()
        self._readRain()
        self._readPaleo()
        self._readFlex()
        self._readOut()
        self._paleoFit()

        self.gravity = 9.81
        self.tNow = self.tStart
        self.saveTime = self.tNow
        if self.strat>0:
            self.saveStrat = self.tNow + self.strat
        else:
            self.saveStrat = self.tEnd + self.tout

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
            self.reflevel = 0

        try:
            meshFile = domainDict['npdata']
        except KeyError as exc:
            print("Key 'npdata' is required and is missing in the 'domain' declaration!")
            raise KeyError('Compressed numpy dataset definition is not defined!')

        if self.reflevel > 0:
            self.meshFile = meshFile+str(self.reflevel)+'.npz'
        else:
            self.meshFile = meshFile+'.npz'

        try:
            with open(self.meshFile) as meshfile:
                pass
        except IOError as exc:
            print("Unable to open numpy dataset: ",self.meshFile)
            raise IOError('The numpy dataset is not found...')

        try:
            self.fast = domainDict['fast']
        except KeyError as exc:
            self.fast = False

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
            self.rStep = timeDict['rstep']
        except KeyError as exc:
            self.rStep = 0

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

        try:
            self.tecStep = timeDict['tec']
        except KeyError as exc:
            self.tecStep = self.tout

        try:
            self.strat = timeDict['strat']
        except KeyError as exc:
            self.strat = 0

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

    def _readTectonic(self):
        """
        Parse tectonic forcing conditions.
        """

        try:
            tecDict = self.input['tectonic']
            tecSort = sorted(tecDict, key=itemgetter('start'))
            for k in range(len(tecSort)):
                tecStart = None
                tMap = None
                zMap = None
                tUniform = None
                tEnd = None
                tStep = None
                try:
                    tecStart = tecSort[k]['start']
                except:
                    print("For each tectonic event a start time is required.")
                    raise ValueError('Tectonic event {} has no parameter start'.format(k))
                try:
                    tMap = tecSort[k]['mapH']+'.npz'
                except:
                    pass
                try:
                    zMap = tecSort[k]['mapV']+'.npz'
                except:
                    pass
                try:
                    tStep = self.tecStep
                except:
                    pass
                try:
                    tEnd = tecSort[k]['end']
                except:
                    pass

                if tMap is not None:
                    if self.meshFile != tMap:
                        try:
                            with open(tMap) as tecfile:
                                pass
                        except IOError as exc:
                            print("Unable to open tectonic file: ",tMap)
                            raise IOError('The tectonic file {} is not found for climatic event {}.'.format(tMap,k))
                else:
                    tMap = 'empty'

                if zMap is not None:
                    if self.meshFile != zMap:
                        try:
                            with open(zMap) as tecfile:
                                pass
                        except IOError as exc:
                            print("Unable to open tectonic file: ",zMap)
                            raise IOError('The tectonic file {} is not found for climatic event {}.'.format(zMap,k))
                else:
                    zMap = 'empty'

                if tMap == 'empty' and zMap == 'empty':
                    print("For each tectonic event a tectonic grid (mapH or mapV) is required.")
                    raise ValueError('Tectonic event {} has no tectonic map (map).'.format(k))

                tmpTec = []
                tmpTec.insert(0, {'start': tecStart, 'tMap': tMap, 'zMap': zMap})

                if k == 0:
                    tecdata = pd.DataFrame(tmpTec, columns=['start', 'tMap', 'zMap'])
                else:
                    tecdata = pd.concat([tecdata,pd.DataFrame(tmpTec, columns=['start', 'tMap', 'zMap'])], ignore_index=True)

                if tStep is not None:
                    if tEnd is not None:
                        tectime = tecStart+tStep
                        while tectime<tEnd:
                            tmpTec = []
                            tmpTec.insert(0, {'start': tectime, 'tMap': tMap, 'zMap': zMap})
                            tecdata = pd.concat([tecdata,pd.DataFrame(tmpTec, columns=['start', 'tMap', 'zMap'])], ignore_index=True)
                            tectime = tectime+tStep

            if tecdata['start'][0] > self.tStart:
                tmpTec = []
                tmpTec.insert(0, {'start': self.tStart, 'tMap': 'empty', 'zMap': 'empty'})
                tecdata = pd.concat([pd.DataFrame(tmpTec, columns=['start', 'tMap', 'zMap']),tecdata], ignore_index=True)
            self.tecdata = tecdata

        except KeyError as exc:
            self.tecdata = None
            pass

        return

    def _readRain(self):
        """
        Parse rain forcing conditions.
        """
        try:
            rainDict = self.input['climate']
            rainSort = sorted(rainDict, key=itemgetter('start'))
            for k in range(len(rainSort)):
                rStart = None
                rUniform = None
                rMap = None
                try:
                    rStart = rainSort[k]['start']
                except:
                    print("For each climate event a start time is required.")
                    raise ValueError('Climate event {} has no parameter start'.format(k))
                try:
                    rUniform = rainSort[k]['uniform']
                except:
                    pass
                try:
                    rMap = rainSort[k]['map']
                except:
                    pass

                if rMap is not None:
                    if self.meshFile != rMap[0]+'.npz':
                        try:
                            with open(rMap[0]+'.npz') as rainfile:
                                pass
                        except IOError as exc:
                            print("Unable to open rain file: ",rMap[0]+'.npz')
                            raise IOError('The rain file {} is not found for climatic event {}.'.format(rMap[0]+'.npz',k))

                        mdata = np.load(rMap[0]+'.npz')
                        rainSet = mdata.files
                    else:
                        mdata = np.load(self.meshFile)
                        rainSet = mdata.files
                    try:
                        rainKey = mdata[rMap[1]]
                    except KeyError as exc:
                        print("Field name {} is missing from rain file {}".format(rMap[1],rMap[0]+'.npz'))
                        print("The following fields are available: ",rainSet)
                        print("Check your rain file fields definition...")
                        raise KeyError('Field name for rainfall is not defined correctly or does not exist!')


                if rMap is None and rUniform is None:
                    print("For each climate event a rainfall value (uniform) or a rainfall grid (map) is required.")
                    raise ValueError('Climate event {} has no rainfall value (uniform) or a rainfall map (map).'.format(k))

                tmpRain = []
                if rMap is None:
                    tmpRain.insert(0, {'start': rStart, 'rUni': rUniform, 'rMap': None, 'rKey': None})
                else:
                    tmpRain.insert(0, {'start': rStart, 'rUni': None, 'rMap': rMap[0]+'.npz', 'rKey': rMap[1]})

                if k == 0:
                    raindata = pd.DataFrame(tmpRain, columns=['start', 'rUni', 'rMap', 'rKey'])
                else:
                    raindata = pd.concat([raindata,pd.DataFrame(tmpRain, columns=['start', 'rUni', 'rMap', 'rKey'])],
                                                         ignore_index=True)

            if raindata['start'][0] > self.tStart:
                tmpRain = []
                tmpRain.insert(0, {'start': self.tStart, 'rUni': 0., 'rMap': None, 'rKey': None})
                raindata = pd.concat([pd.DataFrame(tmpRain, columns=['start', 'rUni', 'rMap', 'rKey']),raindata],
                                                                              ignore_index=True)
            self.raindata = raindata

        except KeyError as exc:
            self.raindata = None
            pass

        return

    def _readPaleo(self):
        """
        Parse paleomap forcing conditions.
        """
        try:
            paleoDict = self.input['paleomap']
            paleoSort = sorted(paleoDict, key=itemgetter('time'))
            for k in range(len(paleoSort)):
                pTime = None
                pMap = None
                try:
                    pTime = paleoSort[k]['time']
                except:
                    print("For each paleomap a given time is required.")
                    raise ValueError('Paleomap {} has no parameter time'.format(k))
                try:
                    pMap = paleoSort[k]['npdata']
                except:
                    pass

                if pMap is not None:

                    try:
                        with open(pMap+'.npz') as meshfile:
                            pass
                    except IOError as exc:
                        print("Unable to open numpy dataset: ",pMap+'.npz')
                        raise IOError('The numpy dataset is not found...')


                tmpPaleo = []
                tmpPaleo.insert(0, {'time': pTime, 'pMap': pMap+'.npz'})

                if k == 0:
                    paleodata = pd.DataFrame(tmpPaleo, columns=['time', 'pMap'])
                else:
                    paleodata = pd.concat([paleodata,pd.DataFrame(tmpPaleo, columns=['time', 'pMap'])],
                                                         ignore_index=True)

            self.paleodata = paleodata
            self.paleoNb = len(paleodata)

        except KeyError as exc:
            self.paleodata = None
            self.paleoNb = 0
            pass

        return

    def _readFlex(self):
        """
        Parse isostatic flexure conditions.
        """
        try:
            flexDict = self.input['flexure']
            self.flexure = True

            try:
                self.dmantle = flexDict['dmantle']
            except KeyError as exc:
                self.dmantle = 3350.0

            try:
                self.dfill = flexDict['dfill']
            except KeyError as exc:
                self.dfill = 2750.0

            try:
                self.young = flexDict['young']
            except KeyError as exc:
                self.young = 1.e11

            try:
                self.fpts = flexDict['npts']
            except KeyError as exc:
                print("The number of points under the influence of a given point load is required.")
                raise ValueError('The nbPts parameter needs to be be defined')

            flexSort = sorted(flexDict, key=itemgetter('time'))
            for k in range(len(flexSort)):
                fTime = None
                fMap = None
                fKey = None
                try:
                    fTime = flexSort[k]['time']
                except:
                    print("For each elastic thickness map a given time is required.")
                    raise ValueError('elastic thickness {} has no parameter time'.format(k))
                try:
                    fMap = flexSort[k]['npdata']
                except:
                    pass

                if fMap is not None:
                    try:
                        with open(fMap+'.npz') as meshfile:
                            pass
                    except IOError as exc:
                        print("Unable to open numpy dataset: ",fMap+'.npz')
                        raise IOError('The numpy dataset is not found...')

                try:
                    fKey = flexSort[k]['key']
                except:
                    pass

                tmpFlex = []
                tmpFlex.insert(0, {'time': fTime, 'fMap': fMap+'.npz', 'fKey': fKey})

                if k == 0:
                    flexdata = pd.DataFrame(tmpFlex, columns=['time', 'fMap', 'fKey'])
                else:
                    flexdata = pd.concat([flexdata,pd.DataFrame(tmpFlex, columns=['time', 'fMap', 'fKey'])],
                                                         ignore_index=True)
            self.flexdata = flexdata

        except KeyError as exc:
            self.flexure = False
            pass

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

        if self.rStep>0:
            self.makedir = False

        return

    def _paleoFit(self):
        """
        Parse paleotopography convergence parameters.
        """

        try:
            paleoDict = self.input['paleofit']
            try:
                self.paleostep = paleoDict['step']
            except KeyError as exc:
                self.paleostep = 10
            try:
                self.cvglimit = paleoDict['cvglimit']
            except KeyError as exc:
                self.paleoDict = 5.

            try:
                self.paleoTopo = paleoDict['paleotopo']
                self.paleoTopo = self.paleoTopo +'.npz'
            except KeyError as exc:
                print("Key 'paleotopo' is required and is missing in the 'paleofit' declaration!")
                raise KeyError('NPZ paleotopography file is not defined!')
            try:
                with open(self.paleoTopo) as paleonpz:
                    pass
            except IOError as exc:
                print("Unable to open NPZ paleo file: ",self.paleoTopo)
                raise IOError('The NPZ paleotopography file is not found...')

            try:
                self.lonlat = paleoDict['lonlat']
                self.lonlat = self.lonlat +'.npz'
            except KeyError as exc:
                print("Key 'lonlat' is required and is missing in the 'paleofit' declaration!")
                raise KeyError('Longitude/latitude coordinates file is not defined!')
            try:
                with open(self.paleoTopo) as paleonpz:
                    pass
            except IOError as exc:
                print("Unable to open lonlat file: ",self.lonlat)
                raise IOError('The lonlat paleotopography file is not found...')

            try:
                self.paleoDir = paleoDict['vdispdir']
            except KeyError as exc:
                print("Key 'vdispdir' is required and is missing in the 'paleofit' declaration!")
                raise KeyError('Paleo-tectonic directory to store vertical displacement is not defined!')

            try:
                os.makedirs(self.paleoDir, exist_ok=True)
            except FileExistsError:
                # Directory already exists
                pass

        except KeyError as exc:
            self.paleostep = 0
            self.cvglimit = 0.
            self.paleoTopo = None
            self.paleoDir = None

        return
