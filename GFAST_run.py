#!/usr/bin/env python
import math
import numpy
import csv
import obspy
import time
import sys
import calendar
import plotting
from datetime import datetime, timedelta
import coord_tools
from eew_data_engine import data_engine_pgd, data_engine_cmt, data_engine_ff, offset_estimator, pgd_estimator
import os
import logging
#import dm_message_writer
#############################################
#GFAST_run.py
#Written by Brendan Crowell, The Ohio State University, August 15, 2025
#This version of the code is designed to build a data buffer from an influxDB
#that stores real-time GNSS positions from Instituto Geofisico. When running the code, user inputs will 
#define the location and timing of the event to process. The station list and location
#is built directly from the available stations within the database at the specific time.
#Additional modules can be added within eew_data_engine
#############################################
if not os.path.exists('output'): #if output folder doesn't exist, make it
        os.makedirs('output')

logger = logging.getLogger("gfast_run")
logger.setLevel(logging.INFO)
if not logger.handlers:
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

logger.info('CLI arguments: %s', sys.argv[:])
lat = sys.argv[1]
lon = sys.argv[2]
dep = sys.argv[3]
timestamp = sys.argv[4]
eqname = sys.argv[5]
ndata = sys.argv[6]
style = sys.argv[7] #make this zero if using the real time server, 1 for playback of events manually added
log_file = os.path.join('output', f'gfast_{eqname}_run.log')
if not any(isinstance(h, logging.FileHandler) for h in logger.handlers):
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
logger.info('Logging to %s', log_file)
logger.info('Selected style=%s for processing', style)
if (int(style) == 0):
        import buffer_init_influxDB

if (int(style) == 1):
	chanfile = sys.argv[8]
	import buffer_init_influxDB_archive
#############################################
#User Input, start with location
#############################################
#print ("Enter latitude:")
#lat = input()
eq_lat = float(lat)
#print ("Enter longitude:")
#lon = input()
eq_lon = float(lon)
#print ("Enter Depth (km, positive down):")
#dep = input()
eq_dep = float(dep)

##########
#Timing
##########
#print ("Enter Origin Time (format, in UTC, is YYYY-MN-DY HR:MN:SC):")
#print ('Enter time as YYYY-MN-DY HR:MN:SC')
#timestamp = input()
struct_time = time.strptime(timestamp, "%Y-%m-%d %H:%M:%S")
ot = datetime.fromtimestamp(time.mktime(struct_time))
logger.info("Event origin time (in UTC): %s", str(ot))
logger.info("Event origin time (unix): %s", str(calendar.timegm(struct_time)))
unixot = int(calendar.timegm(struct_time))
gpstimeoff = (numpy.datetime64('1980-01-06T00:00:00')-numpy.datetime64('1970-01-01T00:00:00'))/ numpy.timedelta64(1, 's')
gpstime = unixot-gpstimeoff
leapsec = coord_tools.gpsleapsec(gpstime)
#Specific to this influxDB is a 5 hour difference between machine time and UTC. I also make the leap second correction
if (int(style) == 0):
        unixotcorr = unixot-18000+leapsec
        logger.info('Using real-time buffer initialization with 5h offset correction. Corrected origin: %s', unixotcorr)
else:
        unixotcorr = unixot+leapsec
        logger.info('Using archive buffer initialization. Corrected origin: %s', unixotcorr)
##########
#EQ Name
##########
#print ("Enter name for earthquake (to be used for output file):")
#eqname = input()
logger.info("Event location (lon,lat,dep(km)): %s, %s, %s", lon, lat, dep)

##########
#Number of seconds to process
##########
#print ("Post earthquake time to output (in seconds):")
#ndata = input()
ndata = int(ndata)
logger.info('Processing window requested: %s seconds', ndata)

##########
#Credentials
##########

from paraminit import Properties
props=Properties('gfast.props')

ip = props.getipaddress()
pt = props.getport()
uname = props.getusername()
pword = props.getpassword()
database = props.getdatabase()
logger.info('InfluxDB target %s:%s database=%s', ip, pt, database)
#############################################
#Run Parameters
#############################################
datarate = 1 #time between samples
logger.info('Sampling interval set to %s second(s)', datarate)
#File names for output
fnamecmt = 'output/gfast_' + eqname + '_cmt.txt' #output file of pgd and cmt results
fnamepgd = 'output/gfast_' + eqname + '_pgd.txt' #output file of pgd and cmt results
fnamefo = 'output/gfast_' + eqname + '_slipmodel_overview.txt' #output file of slip models
fnamefm = 'output/gfast_' + eqname + '_slipmodel.txt' #output file of slip models
fnameff = 'output/gfast_' + eqname + ' _slipfits.txt' #output file of gps vector fits
fnamepgdv =  'output/gfast_' + eqname + '_pgd_values.txt' #output file of pgd values
#Parameters for finite fault inversion
nstr = 10 #number of along strike fault components
ndip = 5 #number of along dip fault components
fcmt = open(fnamecmt,'w')
fpgd = open(fnamepgd,'w')
fpgdv = open(fnamepgdv,'w')
ffo = open(fnamefo,'w')
fff = open(fnameff,'w')
ffm = open(fnamefm,'w')
#############################################
#Build Data buffers
#############################################
if (int(style) == 0):
        logger.info('Building buffers from real-time InfluxDB data')
        [staname,gpst,sta_lat,sta_lon,sta_alt,nbuff,ebuff,ubuff,tvals] = buffer_init_influxDB.data_fromInfluxDB(unixotcorr,unixot,ndata,ip,pt,uname,pword,database)
else:
        logger.info('Building buffers from archive data file %s', chanfile)
        [staname,gpst,sta_lat,sta_lon,sta_alt,nbuff,ebuff,ubuff,tvals] = buffer_init_influxDB_archive.data_fromInfluxDB(unixotcorr,unixot,ndata,chanfile,ip,pt,uname,pword,database)
tbuff = numpy.c_[0:ndata*datarate:datarate] #time buffer

#############################################
#Run GFAST
#############################################
pgd_ver = 0
ffvers = 0
runcmtff = 0
runtime = 0
while runtime < datarate*ndata+1:
        [PGDvals,HYP,EPI,aa] = pgd_estimator(eq_lat,eq_lon,eq_dep,sta_lat,sta_lon,sta_alt,nbuff,ebuff,ubuff,tbuff,runtime)
        OUTPUT_PGD = data_engine_pgd(PGDvals,HYP,EPI,aa,runtime)
	for i in range(0,len(aa)):
		pgdv = "{0:.4f}".format(float(PGDvals[i,0]))
		hypdist = "{0:.2f}".format(float(HYP[i,0])) 
		fpgdv.write(str(runtime)+' '+str(staname[aa[i]])+' '+hypdist+' '+pgdv+'\n')
	MPGD = OUTPUT_PGD[0] #PGD Magnitudes 
	SIG_PGD = OUTPUT_PGD[1] #PGD Sigma
	LEN_PGD = OUTPUT_PGD[2] #Number of stations used for PGD calculation
	VR_PGD = OUTPUT_PGD[3] #PGD Variance Reduction

        if (LEN_PGD > 3):
                mpgd = "{0:.2f}".format(float(MPGD[0,0]))
                vrpgd = "{0:.2f}".format(float(VR_PGD))
                eqlon = "{0:.2f}".format(float(eq_lon))
                eqlat = "{0:.2f}".format(float(eq_lat))
                sigpgd = "{0:.4f}".format(float(SIG_PGD))
                lenpgd = str(LEN_PGD)
                fpgd.write(str(runtime)+' '+mpgd+' '+sigpgd+' '+vrpgd+' '+lenpgd+'\n')
                logger.info('Runtime %s: PGD M=%s sigma=%s VR=%s stations=%s', runtime, mpgd, sigpgd, vrpgd, lenpgd)

        if ((SIG_PGD <= 0.5) or (runcmtff == 1)):
                runcmtff = 1 #if the PGD uncertainty has ever been under 0.5, run the full stack, otherwise skip
                [N,E,U,aa] = offset_estimator(eq_lat,eq_lon,eq_dep,sta_lat,sta_lon,sta_alt,nbuff,ebuff,ubuff,tbuff,runtime)
                SITES = staname[aa]
                OUTPUT_CMT = data_engine_cmt(eq_lat,eq_lon,eq_dep,sta_lat,sta_lon,sta_alt,N,E,U,tbuff,runtime,aa)
	
		MCMT = OUTPUT_CMT[0] #CMT Magnitudes as function of depth
		S1 = OUTPUT_CMT[1] #The 6 moment tensor components as function of depth
		S2 = OUTPUT_CMT[2]
		S3 = OUTPUT_CMT[3]
		S4 = OUTPUT_CMT[4]
		S5 = OUTPUT_CMT[5]
		STR1 = OUTPUT_CMT[6] #The strike, dip and rakes of the main and auxiliary fault planes as a function of depth
		STR2 = OUTPUT_CMT[7]
		DIP1 = OUTPUT_CMT[8]
		DIP2 = OUTPUT_CMT[9]
		RAK1 = OUTPUT_CMT[10]
		RAK2 = OUTPUT_CMT[11]
		VR_CMT = OUTPUT_CMT[12] #Variance reduction for CMT as function of depth
		LEN_CMT = OUTPUT_CMT[13] #Number of stations used for CMT calculation
		dep_vr_cmt = numpy.argmax(VR_CMT) #Depth of greatest variance reduction, CMT

		vrcmt = "{0:.2f}".format(float(VR_CMT[dep_vr_cmt,0]))
		s1 = "{0:.4e}".format(float(S1[dep_vr_cmt,0]))
		s2 = "{0:.4e}".format(float(S2[dep_vr_cmt,0]))
		s3 = "{0:.4e}".format(float(S3[dep_vr_cmt,0]))
		s4 = "{0:.4e}".format(float(S4[dep_vr_cmt,0]))
		s5 = "{0:.4e}".format(float(S5[dep_vr_cmt,0]))
		str1 = "{0:.2f}".format(float(STR1[dep_vr_cmt,0]))
		str2 = "{0:.2f}".format(float(STR2[dep_vr_cmt,0]))
		dip1 = "{0:.2f}".format(float(DIP1[dep_vr_cmt,0]))
		dip2 = "{0:.2f}".format(float(DIP2[dep_vr_cmt,0]))
		rak1 = "{0:.2f}".format(float(RAK1[dep_vr_cmt,0]))
		rak2 = "{0:.2f}".format(float(RAK2[dep_vr_cmt,0]))
		mcmt = "{0:.2f}".format(float(MCMT[dep_vr_cmt,0]))

		lencmt = str(LEN_CMT)
                fcmt.write(str(runtime)+' '+str(dep_vr_cmt)+' '+vrcmt+' '+mcmt+' '+s1+' '+s2+' '+s3+' '+s4+' '+s5+' '+str1+' '+str2+' '+dip1+' '+dip2+' '+rak1+' '+rak2+' '+lencmt+'\n')
                logger.info('Runtime %s: CMT depth index=%s VR=%s M=%s stations=%s', runtime, dep_vr_cmt, vrcmt, mcmt, lencmt)

                if (LEN_CMT > 3):
			str1 = STR1[dep_vr_cmt]
			dip1 = DIP1[dep_vr_cmt]
			str2 = STR2[dep_vr_cmt]
			dip2 = DIP2[dep_vr_cmt]
			mcmt = MCMT[dep_vr_cmt]


                        OUTPUT_FF=data_engine_ff(eq_lat,eq_lon,dep_vr_cmt,mcmt,str1,str2,dip1,dip2,nstr,ndip,sta_lat,sta_lon,sta_alt,N,E,U,tbuff,runtime,aa)
			SSLIP = OUTPUT_FF[0] #Strike-slip along each fault patch
			DSLIP = OUTPUT_FF[1] #Dip-slip along each fault patch
			MFF = OUTPUT_FF[2] #Finite fault magnitude
			EI = OUTPUT_FF[3] #Input east displacements
			NI = OUTPUT_FF[4] #Input north displacements
			UI = OUTPUT_FF[5] #Input up displacements
			EN = OUTPUT_FF[6] #Modeled east displacements
			NN = OUTPUT_FF[7] #Modeled north displacements
			UN = OUTPUT_FF[8] #Modeled up displacements
			STA_LAT = OUTPUT_FF[9] #Station latitudes used in inversion
			STA_LON = OUTPUT_FF[10] #Station longitudes used in inversion
			FAULT_LAT = OUTPUT_FF[11]
			FAULT_LON = OUTPUT_FF[12]
			FAULT_ALT = OUTPUT_FF[13]
			VR_FF1 = OUTPUT_FF[14] #Variance reduction of finite fault inversion for fault plane 1
			VR_FF2 = OUTPUT_FF[15] #Variance reduction of finite fault inversion for fault plane 2
			FaultPlane = OUTPUT_FF[16] #Preferred fault plane, 1 or 2
			FLAT1 = OUTPUT_FF[17]
			FLON1 = OUTPUT_FF[18]
			FLAT2 = OUTPUT_FF[19]
			FLON2 = OUTPUT_FF[20]
			FLAT3 = OUTPUT_FF[21]
			FLON3 = OUTPUT_FF[22]
			FLAT4 = OUTPUT_FF[23]
			FLON4 = OUTPUT_FF[24]


			FFStatus = OUTPUT_FF[25]
			FLEN = OUTPUT_FF[26]
			FWID = OUTPUT_FF[27]

			FDEP1 = OUTPUT_FF[28]
			FDEP2 = OUTPUT_FF[29]
			FDEP3 = OUTPUT_FF[30]
			FDEP4 = OUTPUT_FF[31] 
	
                        if FFStatus == 1:
                                SLIP = numpy.sqrt(numpy.power(SSLIP,2)+numpy.power(DSLIP,2))
				MAXSLIP = numpy.amax(SLIP)
				MAXSLIPloc = numpy.argmax(SLIP)
				MAXSLIPlon = FAULT_LON[MAXSLIPloc,0]
				MAXSLIPlat = FAULT_LAT[MAXSLIPloc,0]
				MAXSLIPdep = FAULT_ALT[MAXSLIPloc,0]
				RAKE = numpy.arctan2(DSLIP,SSLIP)
				RAKEWS = math.degrees(numpy.sum(numpy.multiply(SLIP,RAKE))/numpy.sum(SLIP))
				maxslip = "{0:.4f}".format(float(MAXSLIP))
				maxlon = "{0:.4f}".format(float(MAXSLIPlon))
				maxlat = "{0:.4f}".format(float(MAXSLIPlat))
				maxdep = "{0:.4f}".format(float(MAXSLIPdep))
				rakews = "{0:.4f}".format(float(RAKEWS))
				mff = "{0:.2f}".format(float(MFF))
				vrff1 = "{0:.4f}".format(float(VR_FF1))
				vrff2 = "{0:.4f}".format(float(VR_FF2))
				str11 = "{0:.2f}".format(float(STR1[dep_vr_cmt,0]))
				str22 = "{0:.2f}".format(float(STR2[dep_vr_cmt,0]))
                                ffo.write(str(runtime)+' '+mff+' '+maxlon+' '+maxlat+' '+maxdep+' '+maxslip+' '+rakews+' '+vrff1+' '+vrff2+' '+str11+' '+str22+'\n')
                                logger.info('Runtime %s: FF magnitude=%s max slip=%s km at (%s,%s,%s)', runtime, mff, maxslip, maxlon, maxlat, maxdep)


				l1 = len(SSLIP)
				for i in range (0,l1):
					sslip = "{0:.4f}".format(float(SSLIP[i,0]))
					dslip = "{0:.4f}".format(float(DSLIP[i,0]))
					flat = "{0:.4f}".format(float(FAULT_LAT[i,0]))
					flon = "{0:.4f}".format(float(FAULT_LON[i,0]))
					falt = "{0:.4f}".format(float(FAULT_ALT[i,0]))
					mff = "{0:.4f}".format(float(MFF))
					vrff1 = "{0:.4f}".format(float(VR_FF1))
					vrff2 = "{0:.4f}".format(float(VR_FF2))
					flat1 = "{0:.4f}".format(float(FLAT1[i,0]))
					flon1 = "{0:.4f}".format(float(FLON1[i,0]))
					flat2 = "{0:.4f}".format(float(FLAT2[i,0]))
					flon2 = "{0:.4f}".format(float(FLON2[i,0]))
					flat3 = "{0:.4f}".format(float(FLAT3[i,0]))
					flon3 = "{0:.4f}".format(float(FLON3[i,0]))
					flat4 = "{0:.4f}".format(float(FLAT4[i,0]))
					flon4 = "{0:.4f}".format(float(FLON4[i,0]))
					ffm.write(str(runtime)+' '+str(FaultPlane)+' '+flon+' '+flat+' '+falt+' '+sslip+' '+dslip+' '+mff+' '+vrff1+' '+vrff2+' '+flat1+' '+flon1+' '+flat2+' '+flon2+' '+flat3+' '+flon3+' '+flat4+' '+flon4+'\n')

				l1 = len(E)
				for i in range (0,l1):
					e = "{0:.4f}".format(float(EI[i,0]))
					n = "{0:.4f}".format(float(NI[i,0]))
					u = "{0:.4f}".format(float(UI[i,0]))
					en = "{0:.4f}".format(float(EN[i,0]))
					nn = "{0:.4f}".format(float(NN[i,0]))
					un = "{0:.4f}".format(float(UN[i,0]))
					slat = "{0:.4f}".format(float(STA_LAT[i,0]))
					slon = "{0:.4f}".format(float(STA_LON[i,0]))
					fff.write(str(runtime)+' '+str(SITES[i])+' '+slon+' '+slat+' '+e+' '+n+' '+u+' '+en+' '+nn+' '+un+'\n')

	runtime = runtime+1


if (runcmtff == 0):
        logger.warning('PGD magnitude uncertainty too high for entire run, most likely noise source. No CMT or FF run')
fpgd.close()
fpgdv.close()
fcmt.close()
fff.close()
ffm.close()
ffo.close()

logger.info('Processing completed. Output files stored in output/')

plotting.plotpgdmag(fnamepgd,eqname)
plotting.plotpgdvalues(fnamepgdv,eqname,60)
