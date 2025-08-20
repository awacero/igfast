#!/usr/bin/python
import math
import numpy
import time
import geopy.distance
from coord_tools import ll2utm, utm2ll
from scaling import PGD
from cmt import moment_tensor
from fault_plane import fault_CMT
from RTOkada import rtokada
import matplotlib.pyplot as plt

def offset_estimator(eq_lat,eq_lon,eq_dep,sta_lat,sta_lon,sta_alt,nbuff,ebuff,ubuff,tbuff,runtime):
	coord1 = (float(eq_lat),float(eq_lon))
	epilist = list()
	hyplist = list()
	for i in range(0,len(sta_lon)):
		coord2 = (float(sta_lat[i]),float(sta_lon[i]))
		epi = geopy.distance.distance(coord1,coord2).km
		hyp = math.sqrt(math.pow(epi,2)+math.pow(eq_dep,2))
		hyplist.append(hyp)
		epilist.append(epi)
	epidist = numpy.asarray(epilist).reshape((len(sta_lat),1))
	hypodist = numpy.asarray(hyplist).reshape((len(sta_lat),1))

	effhypodist = hypodist-runtime*0.5
	a1 = numpy.where(effhypodist < (runtime)*2.0)[0]
	if len(a1) > 1:
		N = nbuff[a1,:]
		E = ebuff[a1,:]
		U = ubuff[a1,:]
	else:
		N=0.0
		E=0.0
		U=0.0
	return(N,E,U,a1)

def pgd_estimator(eq_lat,eq_lon,eq_dep,sta_lat,sta_lon,sta_alt,nbuff,ebuff,ubuff,tbuff,runtime):
	coord1 = (float(eq_lat),float(eq_lon))
	epilist = list()
	hyplist = list()
	disp = numpy.sqrt(numpy.power(nbuff,2)+numpy.power(ebuff,2)+numpy.power(ubuff,2))
	for i in range(0,len(sta_lon)):
		coord2 = (float(sta_lat[i]),float(sta_lon[i]))
		epi = geopy.distance.distance(coord1,coord2).km
		hyp = math.sqrt(math.pow(epi,2)+math.pow(eq_dep,2))
		hyplist.append(hyp)
		epilist.append(epi)
	epidist = numpy.asarray(epilist).reshape((len(sta_lat),1))
	hypodist = numpy.asarray(hyplist).reshape((len(sta_lat),1))
	a1 = numpy.where(hypodist < (runtime)*3.0)[0]
	maxD = 0
	if (len(a1) > 0):
		Dnew = disp[a1,0:runtime]
		maxD = numpy.zeros([len(a1),1])
		maxD = numpy.nanmax(Dnew,axis=1,out=maxD,keepdims=True)
	return(maxD,hypodist[a1],epidist[a1],a1)		

def data_engine_pgd(maxD,hypodist,epidist,a1,runtime):
	file99 = numpy.loadtxt('M99.txt')
	t99 = file99[:,0]
	m99 = file99[:,7]
	if len(a1) > 3:
		[MPGD,VR_PGD,IQR] = PGD(100*maxD,hypodist,epidist)
		mpgd = MPGD
		mpgdvr = VR_PGD
		a2 = numpy.where(t99 <= runtime)[0]
		ind99 = a2[len(a2)-1]
		sig = 0.5*math.exp(m99[ind99]-mpgd)
	else:
		mpgd = 0
		mpgdvr = 0
		sig = 99.9
	return(mpgd,sig,len(a1),mpgdvr)



def data_engine_cmt(eq_lat,eq_lon,eq_dep,sta_lat,sta_lon,sta_alt,N,E,U,tbuff,runtime,a1):
        coord1 = (float(eq_lat),float(eq_lon))
        epilist = list()
        hyplist = list()
        for i in range(0,len(sta_lon)):
                coord2 = (float(sta_lat[i]),float(sta_lon[i]))
                epi = geopy.distance.distance(coord1,coord2).km
                hyp = math.sqrt(math.pow(epi,2)+math.pow(eq_dep,2))
                hyplist.append(hyp)
                epilist.append(epi)
        epidist = numpy.asarray(epilist).reshape((len(sta_lat),1))
        hypodist = numpy.asarray(hyplist).reshape((len(sta_lat),1))
        effhypodist = hypodist        
        if len(a1) > 3:
                VARRED = numpy.zeros([50,1])
                L2NRM = numpy.zeros([50,1])
                S1 = numpy.zeros([50,1])
                S2 = numpy.zeros([50,1])
                S3 = numpy.zeros([50,1])
                S4 = numpy.zeros([50,1])
                S5 = numpy.zeros([50,1])
                STR1 = numpy.zeros([50,1])
                STR2 = numpy.zeros([50,1])
                DIP1 = numpy.zeros([50,1])
                DIP2 = numpy.zeros([50,1])
                RAK1 = numpy.zeros([50,1])
                RAK2 = numpy.zeros([50,1])
                MW = numpy.zeros([50,1])
                PERDC = numpy.zeros([50,1])
                VRn = numpy.zeros([50,1])
                VRe = numpy.zeros([50,1])
                VRu = numpy.zeros([50,1])

                ###########################################################################	
                # Grid Search to 50 km. 
                ###########################################################################


                dep = 1
                while dep < 51:
                        [S,VR,Mw,NP1,NP2,PDC] = moment_tensor(sta_lat[a1],sta_lon[a1],sta_alt[a1],N,E,U,eq_lat,eq_lon,dep,runtime,effhypodist[a1],epidist[a1])

                        PERDC[dep-1,0] = PDC*100
                        VARRED[dep-1,0] = VR
                        S1[dep-1,0] = S[0]
                        S2[dep-1,0] = S[1]
                        S3[dep-1,0] = S[2]
                        S4[dep-1,0] = S[3]
                        S5[dep-1,0] = S[4]

                        STR1[dep-1,0] = NP1['strike']
                        STR2[dep-1,0] = NP2['strike']
                        DIP1[dep-1,0] = NP1['dip']
                        DIP2[dep-1,0] = NP2['dip']
                        RAK1[dep-1,0] = NP1['rake']
                        RAK2[dep-1,0] = NP2['rake']

                        MW[dep-1,0] = Mw


                        dep = dep+1

        else:
                VARRED = numpy.zeros([50,1])
                S1 = numpy.zeros([50,1])
                S2 = numpy.zeros([50,1])
                S3 = numpy.zeros([50,1])
                S4 = numpy.zeros([50,1])
                S5 = numpy.zeros([50,1])
                STR1 = numpy.zeros([50,1])
                STR2 = numpy.zeros([50,1])
                DIP1 = numpy.zeros([50,1])
                DIP2 = numpy.zeros([50,1])
                RAK1 = numpy.zeros([50,1])
                RAK2 = numpy.zeros([50,1])
                MW = numpy.zeros([50,1])
                PERDC = numpy.zeros([50,1])
        return(MW,S1,S2,S3,S4,S5,STR1,STR2,DIP1,DIP2,RAK1,RAK2,VARRED,len(a1),PERDC)



def data_engine_ff(eq_lat,eq_lon,eq_dep,MCMT,STR1,STR2,DIP1,DIP2,nstr,ndip,sta_lat,sta_lon,sta_alt,N,E,U,tbuff,runtime,a1):
        coord1 = (float(eq_lat),float(eq_lon))
        epilist = list()
        hyplist = list()
        for i in range(0,len(sta_lon)):
                coord2 = (float(sta_lat[i]),float(sta_lon[i]))
                epi = geopy.distance.distance(coord1,coord2).km
                hyp = math.sqrt(math.pow(epi,2)+math.pow(eq_dep,2))
                hyplist.append(hyp)
                epilist.append(epi)
        epidist = numpy.asarray(epilist).reshape((len(sta_lat),1))
        hypodist = numpy.asarray(hyplist).reshape((len(sta_lat),1))
        effhypodist = hypodist

        if len(a1) > 3:
                [fault_lon1,fault_lat1,fault_alt1,strike1,dip1,dl1,dw1,lon11,lat11,lon12,lat12,lon13,lat13,lon14,lat14,dep11,dep12,dep13,dep14]=fault_CMT(eq_lat,eq_lon,eq_dep,MCMT,STR1,DIP1,nstr,ndip)
                [SSLIP1,DSLIP1,MW1,EN1,NN1,UN1,VR1,Einp,Ninp,Uinp]=rtokada(sta_lat[a1], sta_lon[a1], sta_alt[a1], N, E, U, fault_lat1, fault_lon1, fault_alt1, strike1, dip1, dl1, dw1, nstr, ndip,runtime,effhypodist[a1],eq_lon,eq_lat,MCMT-0.3)

                [fault_lon2,fault_lat2,fault_alt2,strike2,dip2,dl2,dw2,lon21,lat21,lon22,lat22,lon23,lat23,lon24,lat24,dep21,dep22,dep23,dep24]=fault_CMT(eq_lat,eq_lon,eq_dep,MCMT,STR2,DIP2,nstr,ndip)
                [SSLIP2,DSLIP2,MW2,EN2,NN2,UN2,VR2,Einp2,Ninp2,Uinp2]=rtokada(sta_lat[a1], sta_lon[a1], sta_alt[a1], N, E, U, fault_lat2, fault_lon2, fault_alt2, strike2, dip2, dl2, dw2, nstr, ndip,runtime,effhypodist[a1],eq_lon,eq_lat,MCMT-0.3)


                FFStatus = 1
                if (VR1 > VR2):
                        FaultPlane = 1
                        SSLIP = SSLIP1
                        DSLIP = DSLIP1
                        MFF = MW1
                        EN = EN1
                        NN = NN1
                        UN = UN1
                        VR = VR1
                        STR = STR1
                        DIP = DIP1
                        FAULT_LAT = fault_lat1
                        FAULT_LON = fault_lon1
                        FAULT_ALT = fault_alt1
                        FLAT1 = lat11
                        FLON1 = lon11
                        FLAT2 = lat12
                        FLON2 = lon12
                        FLAT3 = lat13
                        FLON3 = lon13
                        FLAT4 = lat14
                        FLON4 = lon14
                        FDEP1 = dep11
                        FDEP2 = dep12
                        FDEP3 = dep13
                        FDEP4 = dep14
                        LEN = dl1
                        WID = dw1
                else:
                        FaultPlane = 2
                        SSLIP = SSLIP2
                        DSLIP = DSLIP2
                        MFF = MW2
                        EN = EN2
                        NN = NN2
                        UN = UN2
                        VR = VR2
                        STR = STR2
                        DIP = DIP2
                        FAULT_LAT = fault_lat2
                        FAULT_LON = fault_lon2
                        FAULT_ALT = fault_alt2
                        FLAT1 = lat21
                        FLON1 = lon21
                        FLAT2 = lat22
                        FLON2 = lon22
                        FLAT3 = lat23
                        FLON3 = lon23
                        FLAT4 = lat24
                        FLON4 = lon24
                        FDEP1 = dep21
                        FDEP2 = dep22
                        FDEP3 = dep23
                        FDEP4 = dep24
                        LEN = dl2
                        WID = dw2
        else:
                SSLIP = -99999
                DSLIP = -99999
                MFF = -99999
                EN = -99999
                NN = -99999
                UN = -99999
                VR = -99999
                STR = -99999
                DIP = -99999
                FAULT_LAT = -99999
                FAULT_LON = -99999
                FAULT_ALT = -99999
                FLAT1 = -99999
                FLON1 = -99999
                FLAT2 = -99999
                FLON2 = -99999
                FLAT3 = -99999
                FLON3 = -99999
                FLAT4 = -99999
                FLON4 = -99999
                FaultPlane = -99999
                Einp = -99999
                Ninp = -99999
                Uinp = -99999
                VR1 = -99999
                VR2 = -99999
                FFStatus = 0
                LEN = -99999
                WID = -99999
                FDEP1=-99999
                FDEP2=-99999
                FDEP3=-99999
                FDEP4=-99999

                
                

        return(SSLIP,DSLIP,MFF,Einp,Ninp,Uinp,EN,NN,UN,sta_lat[a1],sta_lon[a1],FAULT_LAT,FAULT_LON,FAULT_ALT,VR1,VR2,FaultPlane,FLAT1,FLON1,FLAT2,FLON2,FLAT3,FLON3,FLAT4,FLON4,FFStatus,LEN,WID,FDEP1,FDEP2,FDEP3,FDEP4)



