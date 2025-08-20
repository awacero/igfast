from influxdb import InfluxDBClient
import numpy
import time
import calendar
import coord_tools
from datetime import datetime, timedelta
####################################################
####################################################


def data_fromInfluxDB(corrunixot,trueunixot,sec2proc,ip,pt,uname,pword,database):
	client = InfluxDBClient(host=str(ip),port=str(pt), username=str(uname), password=str(pword))
	client.switch_database(str(database))
	#client = InfluxDBClient(host=ip,port=pt, username=uname, password=pword)
	#client.switch_database(database)
	#Convert the origin time to the computer time (5 hours different UTC)
	starttime = (corrunixot-120)*1e9 #I added a 2 minute buffer to the start and stop times
	stoptime = (corrunixot+sec2proc+120)*1e9
	query = "SELECT * FROM gps_position WHERE time >= "+str(int(starttime))+" AND time < "+str(int(stoptime))

	result = client.query(query)
	points = list(result.get_points())
	gpst = numpy.asarray([point['gps_datetime'] for point in points]).astype(float)
	x = numpy.asarray([point['position_x'] for point in points]).astype(float)
	y = numpy.asarray([point['position_y'] for point in points]).astype(float)
	z = numpy.asarray([point['position_z'] for point in points]).astype(float)
	n = numpy.asarray([point['position_n'] for point in points]).astype(float)
	e = numpy.asarray([point['position_e'] for point in points]).astype(float)
	u = numpy.asarray([point['position_u'] for point in points]).astype(float)
	
	latlist=list()
	lonlist=list()
	altlist=list()

	for i in range(0,len(x)):
		[lat,lon,alt]=coord_tools.ecef2lla(x[i],y[i],z[i])
		latlist.append(lat)
		lonlist.append(lon)
		altlist.append(alt)
	lats = numpy.asarray(latlist)
	lons = numpy.asarray(lonlist)
	alts = numpy.asarray(altlist)
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
	aremove = numpy.where(uniquesites == 'TRPG')[0]#removing TRPG because of very poor performance. Delete this and next line if improved.
	uniquesites=numpy.delete(uniquesites,aremove)
	sta_length = len(uniquesites)
	ebuff = numpy.nan*numpy.ones([sta_length,sec2proc])
	nbuff = numpy.nan*numpy.ones([sta_length,sec2proc])
	ubuff = numpy.nan*numpy.ones([sta_length,sec2proc])
	latbuff = numpy.nan*numpy.ones([sta_length,1])
	lonbuff = numpy.nan*numpy.ones([sta_length,1])
	altbuff = numpy.nan*numpy.ones([sta_length,1])
	sitelist=list()
	tbuff = numpy.nan*numpy.ones([sta_length,sec2proc])
	for i in range(0,sta_length):
		a1 = numpy.where(uniquesites[i] == site)[0]
		latbuff[i,0] = lats[a1[0]]
		lonbuff[i,0] = lons[a1[0]]
		altbuff[i,0] = alts[a1[0]]
		sitelist.append(uniquesites[i])
	for i in range(0,sta_length):
		a1 = numpy.where((uniquesites[i] == site) & (toff<=0) & (toff >= -10))[0]
		n0 = numpy.nanmedian(n[a1])
		e0 = numpy.nanmedian(e[a1])
		u0 = numpy.nanmedian(u[a1])
		for j in range(0,sec2proc):
			a1 = numpy.where((uniquesites[i] == site) & (j == toff))[0]
			if (len(a1)>0):
				nbuff[i,j]=float(n[a1[0]])-float(n0)
				ebuff[i,j]=float(e[a1[0]])-float(e0)
				ubuff[i,j]=float(u[a1[0]])-float(u0)
				tbuff[i,j]=j	
	SITES = numpy.asarray(sitelist)
	return(SITES,gpst,latbuff,lonbuff,altbuff,nbuff,ebuff,ubuff,tbuff)


#timestamp = "2025-08-12 23:58:36"
#struct_time = time.strptime(timestamp, "%Y-%m-%d %H:%M:%S")
#ot = datetime.fromtimestamp(time.mktime(struct_time))
#unixot = int(calendar.timegm(struct_time))

#gpstimeoff = (numpy.datetime64('1980-01-06T00:00:00')-numpy.datetime64('1970-01-01T00:00:00'))/ numpy.timedelta64(1, 's')
#gpstime = unixot-gpstimeoff

#leapsec = coord_tools.gpsleapsec(gpstime)
#unixotcorr = unixot-18000+leapsec
#[site,gpst,latbuff,lonbuff,altbuff,nbuff,ebuff,ubuff,toff] = data_fromInfluxDB(unixotcorr,unixot,120,'192.168.132.10','8086','realtime','qsqmpt2004','gpsr')
#print(toff)
#print(nbuff)
#print(n)
#print(toff)
#print(site)
#print(gpst/1e9-unixot-toff)
#print(unixot,timestamp)
#a1 = numpy.where(site == 'MLEC')[0]
#print(nbuff[5,:],ebuff[5,:],ubuff[5,:])
