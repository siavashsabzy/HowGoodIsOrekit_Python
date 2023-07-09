# initiating orekit 11.3.2
import orekit
orekit.initVM()
# importing hipparchus classes ...
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
# importing orekit classes ...
from org.orekit.attitudes import NadirPointing
from org.orekit.bodies import CelestialBodyFactory
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
from org.orekit.forces.gravity import ThirdBodyAttraction
from org.orekit.forces.radiation import IsotropicRadiationSingleCoefficient
from org.orekit.forces.radiation import SolarRadiationPressure
from org.orekit.forces.gravity import Relativity, OceanTides, SolidTides
from org.orekit.forces.drag import IsotropicDrag
from org.orekit.forces.drag import DragForce
from org.orekit.frames import FramesFactory, ITRFVersion
from org.orekit.models.earth.atmosphere.data import MarshallSolarActivityFutureEstimation
from org.orekit.models.earth.atmosphere import DTM2000
from org.orekit.models.earth import ReferenceEllipsoid
from org.orekit.orbits import CartesianOrbit
from org.orekit.propagation import SpacecraftState
from org.orekit.propagation.numerical import NumericalPropagator
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import Constants, IERSConventions, PVCoordinates
from orekit.pyhelpers import setup_orekit_curdir, absolutedate_to_datetime
setup_orekit_curdir("C:\Program Files\Java\jdk-19\lib") # path to "orekit-data.zip" inside brackets, non for current directory
# other packages
import numpy as np
import plotly.express as px
import plotly.io as pio
import plotly.graph_objs as go
import pandas as pd
from bs4 import BeautifulSoup
from datetime import datetime
import webbrowser
webbrowser.register('firefox',None,webbrowser.BackgroundBrowser("C:\\Program Files\\Mozilla Firefox\\firefox.exe"))
pio.renderers.default = "firefox"



# Edited on ........
class SentinelTest_(): # renamed for test ....
    def __init__(self):
        super().__init__()
        self.utc = TimeScalesFactory.getUTC()
        self.eciFrame = FramesFactory.getGCRF()
        self.ecefFrame = FramesFactory.getITRF(ITRFVersion.ITRF_2014, IERSConventions.IERS_2010, True)
        self.RAD2DEG = (180 / np.pi)
        self.DEG2RAD = (np.pi / 180)

    def inter():
        pass
    
    def xlmUTC2OrekitAbsoluteDate(self, DateFromXLMFile):
        parsedTime = datetime.strptime(DateFromXLMFile, "UTC=%Y-%m-%dT%H:%M:%S.%f")
        orekirAbsD = AbsoluteDate(int(parsedTime.year),int(parsedTime.month),int(parsedTime.day), 
                                      int(parsedTime.hour),int(parsedTime.minute), 
                                      (parsedTime.second + parsedTime.microsecond * 1.0e-6),self.utc)
        return orekirAbsD

    def addGroundStation(self, groundStation):
        self.groundStationsList.append(groundStation)

    def addSpacecraft(self, spacecraft):
        self.spacecraftsList.append(spacecraft)

    def getRangeMeasurements(self, propagator, absoluteDateTime, elevationMask):
        pass

    def getPVTFromXLM(self, xlmBlock):
        timeStr = xlmBlock.find_all("UTC")[0].contents[0]
        x = xlmBlock.find_all("X")[0].contents[0]
        y = xlmBlock.find_all("Y")[0].contents[0]
        z = xlmBlock.find_all("Z")[0].contents[0]
        vx = xlmBlock.find_all("VX")[0].contents[0]
        vy = xlmBlock.find_all("VY")[0].contents[0]
        vz = xlmBlock.find_all("VZ")[0].contents[0]
        orekitAbsoluteDate = self.xlmUTC2OrekitAbsoluteDate(timeStr)
        Position = Vector3D( float(x), float(y), float(z))
        Velocity = Vector3D( float(vx), float(vy), float(vz))
        PositionVelocity = PVCoordinates(Position, Velocity)
        return orekitAbsoluteDate, PositionVelocity

    def getInitialOrbitFromXLM(self, xlminit):
        AbsTime, initialPV = self.getPVTFromXLM(xlminit)
        initialEcef2Eci = self.ecefFrame.getTransformTo(self.eciFrame,AbsTime)
        initialInertialPV = initialEcef2Eci.transformPVCoordinates(initialPV)
        initialOrbit = CartesianOrbit(initialInertialPV, self.eciFrame,
                                        AbsTime,
                                        Constants.WGS84_EARTH_MU )
        sentinel = SpacecraftState(initialOrbit)
        return sentinel

        


            
    def getPropagator(self, spacecraft, level, minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance):
        ecefFrame = FramesFactory.getITRF(ITRFVersion.ITRF_2014, IERSConventions.IERS_2010, True)
        eciFrame = FramesFactory.getGCRF()
        moon = CelestialBodyFactory.getMoon()
        sun = CelestialBodyFactory.getSun()
        wgs84Ellipsoid = ReferenceEllipsoid.getWgs84(ecefFrame)
        nadirPointing = NadirPointing(eciFrame, wgs84Ellipsoid)

        thisintegrator = DormandPrince853Integrator(minStep,
                                                     maxStep,
                                                       vecAbsoluteTolerance,
                                                         vecRelativeTolerance)
        
        # determine the level of the propagator
        try:
            if (level == "Low") or (level == "low") or (level == "L") or (level == "l"):
                propagatorCase = 1
            elif (level == "Medium") or (level == "medium") or (level == "M") or (level == "m"):
                propagatorCase = 2
            elif (level == "High") or (level == "high") or (level == "H") or (level == "h"):
                propagatorCase = 3
            else:
                propagatorCase = 0
        except:
            propagatorCase = 0

        if (propagatorCase == 3):

            # Earth gravity field with degree 64 and order 64
            gravityProvider = GravityFieldFactory.getConstantNormalizedProvider(64, 64)
            gravityAttractionModel = HolmesFeatherstoneAttractionModel(ecefFrame, gravityProvider)

            # Third body attraction model
            moon_3dbodyattraction = ThirdBodyAttraction(moon)
            sun_3dbodyattraction = ThirdBodyAttraction(sun)

            # Solar radiation pressure
            isotropicRadiationSingleCoeff = IsotropicRadiationSingleCoefficient(0.02, 1.0)
            solarRadiationPressure = SolarRadiationPressure(sun, wgs84Ellipsoid.getEquatorialRadius(),
                                                            isotropicRadiationSingleCoeff)

            # Relativity
            relativity = Relativity(Constants.EIGEN5C_EARTH_MU)


            oceanicTides = OceanTides(FramesFactory.getITRF(IERSConventions.IERS_2010, True),
                                    Constants.WGS84_EARTH_EQUATORIAL_RADIUS,  Constants.WGS84_EARTH_MU,
                                        5, 5, IERSConventions.IERS_2010, 
                                    TimeScalesFactory.getUT1(IERSConventions.IERS_2010, True))
            solidTidess = SolidTides(FramesFactory.getITRF(IERSConventions.IERS_2010, True),
                                    Constants.WGS84_EARTH_EQUATORIAL_RADIUS,  Constants.WGS84_EARTH_MU,
                                    gravityProvider.getTideSystem(),
                                    IERSConventions.IERS_2010, TimeScalesFactory.getUT1(IERSConventions.IERS_2010, True),
                                    [CelestialBodyFactory.getSun(), CelestialBodyFactory.getMoon()])



            # Atmospheric drag
            #from org.orekit.models.earth.atmosphere import NRLMSISE00
            #atmosphere = NRLMSISE00(msafe, sun, wgs84Ellipsoid)
            msafe = MarshallSolarActivityFutureEstimation(
                MarshallSolarActivityFutureEstimation.DEFAULT_SUPPORTED_NAMES,
                MarshallSolarActivityFutureEstimation.StrengthLevel.AVERAGE)
            atmosphere = DTM2000(msafe, sun, wgs84Ellipsoid)
            isotropicDrag = IsotropicDrag(0.02, 2.2)
            dragForce = DragForce(atmosphere, isotropicDrag)


            satPropagator = NumericalPropagator(thisintegrator)
            satPropagator.setInitialState(spacecraft)
            satPropagator.setAttitudeProvider(nadirPointing)
            satPropagator.addForceModel(gravityAttractionModel)
            satPropagator.addForceModel(moon_3dbodyattraction)
            satPropagator.addForceModel(sun_3dbodyattraction)
            satPropagator.addForceModel(solarRadiationPressure)
            satPropagator.addForceModel(relativity)
            satPropagator.addForceModel(dragForce)
            satPropagator.addForceModel(oceanicTides)
            satPropagator.addForceModel(solidTidess)
        return satPropagator





    def ExBuildWalker(self, num_plane, num_sat, F, refSat):
        # allsat = [[0 for i in range(num_sat)] for i in range(num_plane)]
        pass

    def getECEF(self, time, inertialPV):
        return self.eciFrame.getTransformTo(self.ecefFrame,time).transformPVCoordinates(inertialPV)
    

    def getNormErrors(self, fromfile, predicted):
        x_resi = (fromfile.getPosition().getX() - predicted.getPosition().getX())
        y_resi = (fromfile.getPosition().getY() - predicted.getPosition().getY())
        z_resi = (fromfile.getPosition().getZ() - predicted.getPosition().getZ())
        norm_resi = (x_resi**2 + y_resi**2 + z_resi**2)**0.5

        vx_resi = (fromfile.getVelocity().getX() - predicted.getVelocity().getX())
        vy_resi = (fromfile.getVelocity().getY() - predicted.getVelocity().getY())
        vz_resi = (fromfile.getVelocity().getZ() - predicted.getVelocity().getZ())
        vnorm_resi = (vx_resi**2 + vy_resi**2 + vz_resi**2)**0.5

        return norm_resi, vnorm_resi



    def getSentinelError(self, path2file):
        fullDataFrame = pd.DataFrame(columns=['Position [m]',
                                              'Velocity [m/s]', 
                                              'Predicted Position [m]',
                                              'Predicted Velocity [m/s]',
                                              'Position Prediction Error [m]',
                                              'Velocity Prediction Error [m/s]'])
        with open(path2file, 'r') as ephfile:
            sentineldata = ephfile.read()
        Bs_sentineldata = BeautifulSoup(sentineldata, "xml")
        fileBlocks = Bs_sentineldata.find_all('OSV')
        initialEphTime = self.xlmUTC2OrekitAbsoluteDate(fileBlocks[0].find_all("UTC")[0].contents[0])
        endEphTime = self.xlmUTC2OrekitAbsoluteDate(fileBlocks[len(fileBlocks) - 1].find_all("UTC")[0].contents[0])
        initialOrbit = self.getInitialOrbitFromXLM(fileBlocks[0])
        sentinelpropagator = self.getPropagator(initialOrbit, "H", 0.1,1.0,0.01,0.01)
        blockIndex = int(0)
        currentDateTime = initialEphTime
        print("The first block is going to be ignored")
        while (currentDateTime.compareTo(endEphTime) < 0):
            blockIndex += 1
            print(blockIndex," of ", len(fileBlocks)-1)
            currentDateTime = self.xlmUTC2OrekitAbsoluteDate(fileBlocks[blockIndex].find_all("UTC")[0].contents[0])
            currentInertialPV = sentinelpropagator.propagate(currentDateTime).getPVCoordinates()
            currentFixedPV = self.getECEF(currentDateTime,currentInertialPV)
            timeFromFile, fixedPVFromFile = self.getPVTFromXLM(fileBlocks[blockIndex])
            if (currentDateTime.compareTo(timeFromFile) != 0):
                break
            positionNormError, VelocityNormError = self.getNormErrors(fixedPVFromFile, currentFixedPV)
            fullDataFrame.loc[absolutedate_to_datetime(currentDateTime)] = (
                fixedPVFromFile.getPosition(),
                fixedPVFromFile.getVelocity(),
                currentFixedPV.getPosition(),
                currentFixedPV.getVelocity(),
                positionNormError,
                VelocityNormError)
        return fullDataFrame




pathTo_S1A_OPER_AUX_RESORB_OPOD_File = 'S1A_OPER_AUX_RESORB_OPOD_20210612T161045_V20210612T121457_20210612T153227.EOF'



fullDataFrame = SentinelTest_().getSentinelError(pathTo_S1A_OPER_AUX_RESORB_OPOD_File)




trace = go.Scattergl(
    x=fullDataFrame.index, y=fullDataFrame['Position Prediction Error [m]'],
    mode='markers',
    name='Position Prediction Error [m]')
layout = go.Layout(
    title = 'Position residuals',
    xaxis = dict(title = 'Datetime UTC'),
    yaxis = dict(title = 'Position Prediction Error [m]'))
fig = dict(data=trace, layout=layout)
pio.show(fig)