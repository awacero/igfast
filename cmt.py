#!/usr/bin/python
import math
import numpy
import okadapoint
import obspy.imaging.beachball
import obspy.signal
import scipy
from pyproj import Transformer

def to_utm(lat,lon):
	transformer = Transformer.from_crs("EPSG:4326","EPSG:32717")
	utm_x, utm_y = transformer.transform(lat,lon)
	return (utm_x, utm_y)

def moment_tensor(sta_lat, sta_lon, sta_alt, n, e, u, eq_lat, eq_lon, eq_alt,runtime,effhypodist,repi):
        l1=len(sta_lat)
        l2=1

        xrs = numpy.zeros([l1,1])
        yrs = numpy.zeros([l1,1])
        zrs = numpy.zeros([l1,1])
        for j in range (0, l1):
                (x1, y1) = to_utm(float(sta_lat[j]),float(sta_lon[j]))
                (x2, y2) = to_utm(float(eq_lat),float(eq_lon))	
                #(x1,y1) = ll2utm(sta_lon[j],sta_lat[j]+25, eq_lon, eq_lat+25)
                #(x2,y2) = ll2utm(eq_lon,eq_lat+25, eq_lon, eq_lat+25)
                xrs[j,0] = (x1-x2)
                yrs[j,0] = (y1-y2)
                zrs[j,0] = (sta_alt[j]+eq_alt*1000)

        U = numpy.zeros([3*l1,1])

        W = numpy.zeros([3*l1,3*l1])
        for i in range (0, l1):
                #efftime = math.ceil(effhypodist[i]/3)
                NN = n[i,runtime-10:runtime]
                EE = e[i,runtime-10:runtime]
                UU = u[i,runtime-10:runtime]
                N = numpy.nanmean(NN)
                E = numpy.nanmean(EE)
                Z = numpy.nanmean(UU)

                W[3*i,3*i]= 1/1
                W[3*i+1,3*i+1]= 1/1
                W[3*i+2,3*i+2]= 1/5
                        
                U[3*i,0]= N
                U[3*i+1,0]= E
                U[3*i+2,0]= -Z
                nout = "{0:.4f}".format(float(N))
                eout = "{0:.4f}".format(float(E))
                zout = "{0:.4f}".format(float(Z))
                latout = "{0:.4f}".format(float(sta_lat[i]))
                lonout = "{0:.4f}".format(float(sta_lon[i]))


        G = okadapoint.greenF(yrs, xrs, -zrs) #Compute Green's function
        S = numpy.linalg.lstsq(numpy.dot(W,G),numpy.dot(W,U),rcond=None)[0]
        M12 = S[0,0]
        M13 = S[1,0]
        M33 = S[2,0]
        M23 = S[4,0]
        M11 = S[3,0]-0.5*S[2,0]
        M22 = -S[3,0]-0.5*S[2,0]

        MNE = S[0,0]
        MND = S[1,0]
        MDD = S[2,0]
        MED = S[4,0]
        MNN = S[3,0]-0.5*S[2,0]
        MEE = -S[3,0]-0.5*S[2,0]

        M_devi = numpy.array([[M11,M12,M13],[M12,M22,M23],[M13,M23,M33]])
        eigenwtot, eigenvtot = numpy.linalg.eig(M_devi)
        eigenw1, eigenv1 = numpy.linalg.eig(M_devi)
        eigenw = numpy.real(numpy.take(eigenw1, numpy.argsort(abs(eigenwtot))))
        eigenv = numpy.real(numpy.take(eigenv1, numpy.argsort(abs(eigenwtot)), 1))
        eigenw_devi = numpy.real(numpy.take(eigenw1, numpy.argsort(abs(eigenw1))))
        eigenv_devi = numpy.real(numpy.take(eigenv1, numpy.argsort(abs(eigenw1)), 1))
        M0_devi = max(abs(eigenw_devi))
        F = -eigenw_devi[0] / eigenw_devi[2]
        M_DC_percentage = (1 - 2 * abs(F))

        Mo = math.pow(math.pow(M11,2)+math.pow(M22,2)+math.pow(M33,2)+2*math.pow(M12,2)+2*math.pow(M13,2)+2*math.pow(M23,2),0.5)/math.pow(2,0.5)
        if (Mo == 0):
                Mw = 0
        else:
                Mw = 2*math.log10(Mo)/3-6.03

        UP = numpy.dot(G,S)

        dms = U-UP
        Emfit = numpy.zeros([l1,1])
        Nmfit = numpy.zeros([l1,1])
        Umfit = numpy.zeros([l1,1])
        Efor = numpy.zeros([l1,1])
        Nfor = numpy.zeros([l1,1])
        Ufor = numpy.zeros([l1,1])
        for i in range (0, l1):
                Nmfit[i,0] = dms[3*i,0]
                Emfit[i,0] = dms[3*i+1,0]
                Umfit[i,0] = dms[3*i+2,0]
                Nfor[i,0] = U[3*i,0]
                Efor[i,0] = U[3*i+1,0]
                Ufor[i,0] = U[3*i+2,0]
                
        VR = 100*(1-numpy.sum(numpy.sqrt((dms)**2))/numpy.sum(numpy.sqrt((U)**2)))

        mt = obspy.imaging.beachball.MomentTensor(M33,M11,M22,M13,-M23,-M12,26)


        axes = obspy.imaging.beachball.mt2axes(mt)
        plane1 = obspy.imaging.beachball.mt2plane(mt)
        plane2 = obspy.imaging.beachball.aux_plane(plane1.strike,plane1.dip,plane1.rake)
        T = {'azimuth':axes[0].strike,'plunge':axes[0].dip}
        N = {'azimuth':axes[1].strike,'plunge':axes[1].dip}
        P = {'azimuth':axes[2].strike,'plunge':axes[2].dip}
        NP1 = {'strike':plane1.strike,'dip':plane1.dip,'rake':plane1.rake}
        NP2 = {'strike':plane2[0],'dip':plane2[1],'rake':plane2[2]}


        return(S, VR, Mw, NP1, NP2, M_DC_percentage)
