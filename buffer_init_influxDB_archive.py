from influxdb import InfluxDBClient
import numpy
import time
import calendar
import coord_tools
import matplotlib.pylab as plt
from datetime import datetime, timedelta
####################################################
#Credentials
####################################################


def data_fromInfluxDB(corrunixot,trueunixot,sec2proc,chanfile,ip,pt,uname,pword,database):
	A = numpy.loadtxt(chanfile,dtype=str)
	lat_chan = A[:,4].astype(float)
	lon_chan = A[:,5].astype(float)
	site_chan = A[:,1].astype(str)
	client = InfluxDBClient(host=str(ip),port=str(pt), username=str(uname), password=str(pword))
	client.switch_database(str(database))
	#Convert the origin time to the computer time (5 hours different UTC)
	starttime = (corrunixot-120)*1e9 #I added a 2 minute buffer to the start and stop times
	stoptime = (corrunixot+sec2proc+120)*1e9
	query = "SELECT * FROM gps_position WHERE time >= "+str(int(starttime))+" AND time < "+str(int(stoptime))

	result = client.query(query)
	points = list(result.get_points())
	gpst = numpy.asarray([point['gps_datetime'] for point in points]).astype(float)
	n = numpy.asarray([point['position_n'] for point in points]).astype(float)
	e = numpy.asarray([point['position_e'] for point in points]).astype(float)
	u = numpy.asarray([point['position_u'] for point in points]).astype(float)
	

	try:
		site = numpy.asarray([point['site_id'] for point in points])
		uniquesites = numpy.unique(site)
	except Exception:
		try:
			site = numpy.asarray([point['site'] for point in points])
			uniquesites = numpy.unique(site)
		except Exception:
			print("no sites found")
	time = numpy.asarray([point['time'] for point in points])
	toff = gpst/1e9-trueunixot
	uniquesites = numpy.unique(site)
	#print(uniquesites)
	sta_length = len(site_chan)
	ebuff = numpy.nan*numpy.ones([sta_length,sec2proc])
	nbuff = numpy.nan*numpy.ones([sta_length,sec2proc])
	ubuff = numpy.nan*numpy.ones([sta_length,sec2proc])
	latbuff = numpy.nan*numpy.ones([sta_length,1])
	lonbuff = numpy.nan*numpy.ones([sta_length,1])
	altbuff = numpy.nan*numpy.ones([sta_length,1])
	sitelist=list()
	tbuff = numpy.nan*numpy.ones([sta_length,sec2proc])
	for i in range(0,sta_length):
		#a1 = numpy.where(uniquesites[i] == site_chan)[0]
		#print(uniquesites[i],a1)
		latbuff[i,0] = lat_chan[i]
		lonbuff[i,0] = lon_chan[i]
		altbuff[i,0] = 0
		sitelist.append(site_chan[i])
	for i in range(0,sta_length):
		a1 = numpy.where((site_chan[i] == site) & (toff<=0) & (toff >= -10))[0]
		n0 = numpy.nanmedian(n[a1])
		e0 = numpy.nanmedian(e[a1])
		u0 = numpy.nanmedian(u[a1])
		for j in range(0,sec2proc):
			a1 = numpy.where((site_chan[i] == site) & (j == toff))[0]
			if (len(a1)>0):
				nbuff[i,j]=(float(n[a1[0]])-float(n0))/1e6
				ebuff[i,j]=(float(e[a1[0]])-float(e0))/1e6
				ubuff[i,j]=(float(u[a1[0]])-float(u0))/1e6
				tbuff[i,j]=j	
	SITES = numpy.asarray(sitelist)
	return(SITES,gpst,latbuff,lonbuff,altbuff,nbuff,ebuff,ubuff,tbuff)


#timestamp = "2016-04-16 23:58:36"
#struct_time = time.strptime(timestamp, "%Y-%m-%d %H:%M:%S")
#ot = datetime.fromtimestamp(time.mktime(struct_time))
#unixot = int(calendar.timegm(struct_time))

#gpstimeoff = (numpy.datetime64('1980-01-06T00:00:00')-numpy.datetime64('1970-01-01T00:00:00'))/ numpy.timedelta64(1, 's')
#gpstime = unixot-gpstimeoff

#leapsec = coord_tools.gpsleapsec(gpstime)
#unixotcorr = unixot+leapsec


#[site,gpst,latbuff,lonbuff,altbuff,nbuff,ebuff,ubuff,toff] = data_fromInfluxDB(unixotcorr,unixot,180,'Ecuador2016_disp_pgd_v2.chan')
#print(toff)
#print(nbuff)
#print(toff)
#print(site,latbuff,lonbuff)
#print(unixot,timestamp)


#print(site)
#a1 = numpy.where(site == 'IBEC')[0]

#plt.plot(nbuff[2,:])
#plt.plot(ebuff[2,:])
#plt.plot(ubuff[2,:])
#plt.savefig('cabp.png')
