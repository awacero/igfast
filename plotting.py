import matplotlib.pylab as plt
import numpy


def plotpgdmag(fname,eqname):
	a = numpy.loadtxt(fname)
	if (len(a) > 0):
		t = a[:,0].astype(float)
		m = a[:,1].astype(float)
		munc = a[:,2].astype(float)
		vr = a[:,3].astype(float)
		numsta = a[:,4].astype(float)


		a1 = numpy.where(munc < 0.5)[0]
		a2 = numpy.where(munc >= 0.5)[0]
		plt.subplot(3,1,1)
		if (len(a1) > 0):
			plt.plot(t[a1],m[a1],'o',color='green')
		if (len(a2) > 0):
			plt.plot(t[a2],m[a2],'s',color='red')
		plt.ylabel('Magnitude')
		plt.title(str(eqname))
		plt.subplot(3,1,2)
		plt.plot(t,vr)
		plt.ylabel('VR (%)')
		plt.subplot(3,1,3)
		plt.plot(t,numsta)
		plt.xlabel('Time (s)')
		plt.ylabel('# Stations')

		outfile = 'output/'+eqname+'_pgdevolution.png'
		plt.savefig(str(outfile), dpi=450)
		plt.clf()
		plt.close()
	return


def plotpgdvalues(fname,eqname,tplot):	
	file99 = numpy.loadtxt('M99.txt')
	t99 = file99[:,0]
	m99 = file99[:,7]
	a2 = numpy.where(t99 <= tplot)[0]
	ind99 = a2[len(a2)-1]
	m99val = m99[ind99]

	a = numpy.loadtxt(fname,dtype=str)
	if (len(a) > 0):
		t = a[:,0].astype(float)
		site = a[:,1]
		hyp = a[:,2].astype(float)
		pgd = a[:,3].astype(float)
		Mpred = (numpy.log10(pgd*100)+3.841)/(0.937-0.127*numpy.log10(hyp))
		a1 = numpy.where(t == tplot)[0]

		plt.subplot(1,2,1)
		plt.loglog(hyp[a1],pgd[a1]*100,'o')
		plt.xlim((1,500))
		plt.ylim((0.1,500))
		plt.xlabel('Hypocentral Distance (km)')
		plt.ylabel('PGD (cm)')
		plt.title(str(eqname)+', '+str(tplot)+' seconds')
		plt.subplot(1,2,2)
		plt.semilogy(Mpred[a1],pgd[a1]*100,'o')
		plt.vlines(x=[m99val],ymin=0.1,ymax=500,colors=['r'])
		plt.ylim((0.1,500))
		plt.xlabel('Predicted Magnitude')
		outfile = 'output/'+eqname+'_pgdvals_'+str(tplot)+'.png'
		plt.savefig(str(outfile),dpi=450)
		plt.clf()
		plt.close()
	return


#plotpgdvalues('output/gfast_test_pgd_values.txt','test',60)
