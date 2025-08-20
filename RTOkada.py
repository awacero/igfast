#!/usr/bin/python
import math
import numpy
import okadagreen
from pyproj import Transformer

def to_utm(lat,lon):
        transformer = Transformer.from_crs("EPSG:4326","EPSG:32717")
        utm_x, utm_y = transformer.transform(lat,lon)
        return (utm_x, utm_y)

def rtokada(sta_lat, sta_lon, sta_alt, n, e, u, fault_lat, fault_lon, fault_alt, strike, dip, LEN, WID, nstr, ndip, runtime, effhypodist,eq_lon, eq_lat, MCMT):
        l1=len(sta_lat)
        l2=len(fault_lat)

        LEN = LEN*1000 
        WID = WID*1000
        xrs = numpy.zeros([l2,l1])
        yrs = numpy.zeros([l2,l1])
        zrs = numpy.zeros([l2,l1])

        for i in range (0, l2):
                for j in range (0, l1):			
                        (x1,y1) = to_utm(sta_lat[j],sta_lon[j])
                        (x2,y2) = to_utm(fault_lat[i],fault_lon[i])
                        xrs[i,j] = (x1-x2)
                        yrs[i,j] = (y1-y2)
                        zrs[i,j] = sta_alt[j]+fault_alt[i]*1000

        G = okadagreen.greenF(xrs, yrs, zrs, strike, dip, WID, LEN) #Compute Green's functions
        T = numpy.zeros([(2*ndip*nstr),2*l2]) #Prefill matrix with zeros
        TU = numpy.zeros([(2*ndip*nstr),1]) #Appended to observation vector. This minimizes the difference between adjacent slip cells

        k=0
        for j in range (0, ndip):
                for i in range (0, nstr):
                        for m in range (0, 2):
                                index1 = j*nstr+i
                                index2 = j*nstr+i-1
                                index3 = j*nstr+i+1
                                index4 = (j-1)*nstr+i
                                index5 = (j+1)*nstr+i

                                if (index1 >= 0 and index1 < l2):
                                        T[k,2*index1+m] = -2.0/LEN[0]/LEN[0]*1000*1000-2.0/WID[0]/WID[0]*1000*1000
                                if (index2 >= 0 and index2 < l2):
                                        T[k,2*index2+m] = 1.0/LEN[0]/LEN[0]*1000*1000
                                if (index3 >= 0 and index3 < l2):
                                        T[k,2*index3+m] = 1.0/LEN[0]/LEN[0]*1000*1000
                                if (index4 >= 0 and index4 < l2):
                                        T[k,2*index4+m] = 1.0/WID[0]/WID[0]*1000*1000
                                if (index5 >= 0 and index5 < l2):
                                        T[k,2*index5+m] = 1.0/WID[0]/WID[0]*1000*1000

                                k=k+1
        W = numpy.zeros([3*l1,3*l1])
        U = numpy.zeros([3*l1,1]) #Create data vector
        Einp = numpy.zeros([l1,1]) #Create data vector
        Ninp = numpy.zeros([l1,1]) #Create data vector
        Uinp = numpy.zeros([l1,1]) #Create data vector
        for i in range (0, l1):
                NN = n[i,runtime-10:runtime]
                EE = e[i,runtime-10:runtime]
                UU = u[i,runtime-10:runtime]

                Einp[i,0] = numpy.nanmean(EE)
                Ninp[i,0] = numpy.nanmean(NN)
                Uinp[i,0] = numpy.nanmean(UU)

                U[3*i,0]= numpy.nanmean(EE)
                U[3*i+1,0]= numpy.nanmean(NN)
                U[3*i+2,0]= numpy.nanmean(UU)
                W[3*i,3*i]= 1/1
                W[3*i+1,3*i+1]=1/1
                W[3*i+2,3*i+2]=1/5
        #print(U)
        UD = numpy.vstack((numpy.dot(W,U),TU))
        #UD = numpy.vstack((U,TU)) 
        if float(MCMT[0]) >= 8.5:
                SFac = 1/0.15
        if float(MCMT[0]) >= 8.0 and float(MCMT[0]) < 8.5:
                SFac = 1/0.5
        if float(MCMT[0]) >= 7.5 and float(MCMT[0]) < 8.0:
                SFac = 1.0
        if float(MCMT[0]) < 7.5:
                SFac = 0.5
                        
        lampred = 1.0/math.pow(l2*2,2)/numpy.mean(numpy.absolute(G))/20*LEN[0]*WID[0]/1000/1000
        T2 = T*lampred

        G2 = numpy.vstack((numpy.dot(W,G),T2))
        S = numpy.linalg.lstsq(G2,UD,rcond=None)[0]
        #print(S)
        UP = numpy.dot(G,S) #Forward model for model fits
        VR = (1-numpy.linalg.norm(UP-U)**2/numpy.linalg.norm(U)**2)*100 #variance reduction

        #VR = numpy.linalg.norm(UP-U)

        SSLIP = numpy.zeros([l2,1])
        DSLIP = numpy.zeros([l2,1])

        for i in range (0,l2):
                SSLIP[i,0] = S[2*i,0]
                DSLIP[i,0] = S[2*i+1,0]

        EN = numpy.zeros([l1,1])
        NN = numpy.zeros([l1,1])
        UN = numpy.zeros([l1,1])
        for i in range (0, l1):
                EN[i,0]=UP[3*i,0]
                NN[i,0]=UP[3*i+1,0]
                UN[i,0]=UP[3*i+2,0]

        ST = numpy.sqrt(SSLIP**2+DSLIP**2)

        Mo = numpy.sum(3.0e10*ST*LEN*WID)
        if (Mo > 0):
                MW = 2.0/3.0*math.log10(Mo/1.0e-7)-10.7
        else:
                MW = 0

	
        return(SSLIP,DSLIP,MW,EN,NN,UN,VR,Einp,Ninp,Uinp)
