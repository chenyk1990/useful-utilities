## This is a script for the test of deblending using sparsity promotion with shaping regularization.

## Copyright (C) 2013
## Yangkang Chen, The Univeristy of Texas at Austin

## The are several assumptions for implementing the following scripts:
## 1. data1 and data2 have the same dimensions 
## 2. data1 and data2 are sparse in fourier or seislet domains
## 3. the unblended data is known (in other words it's used for testing). 
##    Otherwise, the unblended1(2) is substituted by blended1(2), in which case
##    The error diagram has no physical meanings.
## 4. there is no self interference

from rsf.proj import*
def Grey(data,other): 
	Result(data,'grey label2=Trace unit2="" label1=Time unit1="s" title="" wherexlabel=b wheretitle=t color=b %s'%other)

def Greyplot(data,other): 
	Plot(data,'grey label2=Trace unit2="" label1=Time unit1="s" title="" wherexlabel=b wheretitle=t color=b %s'%other)

def Graph(data,other):
	Result(data,'graph label1="" label2="" unit1="" unit2=""  title="" wherexlabel=b wheretitle=t %s' %other)

def Graphplot(data,other):
	Plot(data,'graph label1="" label2="" unit1="" unit2=""  title="" wherexlabel=b wheretitle=t %s' %other)

def shapeslet(sig1,sig2,dip1,dip2,mute1,mute2,padno,n2,mode,thr,i): #mute1(2) is a string, when there is no mute, mute1(2) is "cp"
    Flow(sig1+"p%d"%(i+1),[sig1,dip1],
	 '''
	 %s |pad n2=%d  |
         seislet dip=${SOURCES[1]} eps=0.1 
		  adj=y inv=y unit=y type=b |
         threshold1 type=%s ifperc=1 thr=%f |
         seislet dip=${SOURCES[1]} eps=0.1  
			inv=y unit=y type=b | 
         window n2=%d |%s
         '''%(mute1,padno,mode,thr,n2,mute1))
    Flow(sig2+"p%d"%(i+1),[sig2,dip2],
	 ''' 
	 %s|pad n2=%d  |
         seislet dip=${SOURCES[1]} eps=0.1 
		  adj=y inv=y unit=y type=b |
         threshold1 type=%s ifperc=1 thr=%f |
         seislet dip=${SOURCES[1]} eps=0.1  
			inv=y unit=y type=b | 
         window n2=%d | %s
         '''%(mute2,padno,mode,thr,n2,mute2))
    return (sig1+"p%d"%(i+1),sig2+"p%d"%(i+1))

def shapefft(sig1,sig2,mute1,mute2,mode,thr,i):#mute1(2) is a string, when there is no mute, mute1(2) is "cp"
    Flow(sig1+"p%d"%(i+1),sig1,
	 '''
	 %s|fft1 | fft3 axis=2 pad=1|
         threshold1 type=%s ifperc=1 thr=%f |
         fft3 inv=y axis=2 | fft1 inv=y |%s
         '''%(mute1,mode,thr,mute1))
    Flow(sig2+"p%d"%(i+1),sig2,
	 '''
	 %s|fft1 | fft3 axis=2 pad=1|
         threshold1 type=%s ifperc=1 thr=%f |
         fft3 inv=y axis=2 | fft1 inv=y | %s
         '''%(mute2,mode,thr,mute2))
    return (sig1+"p%d"%(i+1),sig2+"p%d"%(i+1))

def step(sig1,sig2,blended1,blended2,init1,init2,shottime1,shottime2,i,fraction):
    Flow(init1+"%d"%(i+1),[sig2,shottime1,shottime2,blended1,sig1], 
         '''
          blend shot_time_in=${SOURCES[2]} shot_time_out=${SOURCES[1]}
	    | add scale=%f,%f,%f ${SOURCES[3]} ${SOURCES[4]}
         '''%(-fraction,fraction,1-fraction))

    Flow(init2+"%d"%(i+1),[sig1,shottime2,shottime1,blended2,sig2],
         '''
          blend shot_time_in=${SOURCES[2]} shot_time_out=${SOURCES[1]}
	    | add scale=%f,%f,%f ${SOURCES[3]} ${SOURCES[4]}
         '''% (-fraction,fraction,1-fraction))
    return (init1+"%d"%(i+1),init2+"%d"%(i+1))

def deblendfft(unblended1,
	  unblended2,
	  blended1,
          blended2,
	  init1,
	  init2,
	  deblended1,
	  deblended2,
	  shottime1,
	  shottime2,
	  mute1,
	  mute2,
	  n1,
	  n2,
	  niter,
	  mode,
	  thr,
	  pc,
	  fraction):

    sig1, sig2 = (init1,init2)
    sigs1=[]
    sigs2=[]
    diffsa=[]
    diffsb=[]
    diffsc=[]
    diffsd=[]
    for i in range(niter):
 	old1,old2 = (sig1,sig2)
	diffa=init1+'-diffa%d'%i
	diffb=init2+'-diffb%d'%i
	diffc=init1+'-diffc%d'%i
	diffd=init2+'-diffd%d'%i

	(nsig1,nsig2)=step(sig1,sig2,blended1,blended2,init1,init2,shottime1,shottime2,i,fraction)
	(sig1,sig2)=shapefft(nsig1,nsig2,mute1,mute2,mode,thr,i)

	Flow(diffa,[old1,sig1],'diff match=${SOURCES[1]} | math output="input/%d/%d"'%(n1,n2))
	Flow(diffb,[old2,sig2],'diff match=${SOURCES[1]} | math output="input/%d/%d"'%(n1,n2))
        Flow(diffc,[unblended1,sig1],'diff match=${SOURCES[1]} | math output="input/%d/%d"'%(n1,n2))
	Flow(diffd,[unblended2,sig2],'diff match=${SOURCES[1]} | math output="input/%d/%d"'%(n1,n2))
	
    	Greyplot(sig1,'title="Esti R 1 (iter=%d)" pclip=%d'% (i+1,pc))
    	Greyplot(sig2,'title="Esti R 2 (iter=%d)" pclip=%d' % (i+1,pc))	

	sigs1.append(sig1)
    	sigs2.append(sig2)
        diffsa.append(diffa)
        diffsb.append(diffb)
	diffsc.append(diffc)
	diffsd.append(diffd)

    Flow(init1+'-diffsa',diffsa,'cat axis=2 ${SOURCES[1:%d]}|math output="log(input)/log(10)" | transp | put n1=%d o1=1 d1=1 '%(len(diffsa),len(diffsa)))
    Flow(init2+'-diffsb',diffsb,'cat axis=2 ${SOURCES[1:%d]}|math output="log(input)/log(10)" | transp |  put n1=%d o1=1 d1=1'%(len(diffsb),len(diffsb)))
    Flow(init1+'-diffsc',diffsc,'cat axis=2 ${SOURCES[1:%d]}|math output="log(input)/log(10)" | transp |  put n1=%d o1=1 d1=1'%(len(diffsc),len(diffsc)))
    Flow(init2+'-diffsd',diffsd,'cat axis=2 ${SOURCES[1:%d]}|math output="log(input)/log(10)" | transp |  put n1=%d o1=1 d1=1'%(len(diffsd),len(diffsd)))

    Graph(init1+'-diffsa','title="Data error 1" label1=Iteration label2="c(n)" symbosz=5 symbol="o" min1=0 max1=%d d1=1'%niter)	
    Graph(init2+'-diffsb','title="Data error 2" label1=Iteration label2="c(n)" symbosz=5 symbol="o" min1=0 max1=%d d1=1'%niter)	
    Graph(init1+'-diffsc','title="Model error 1" label1=Iteration label2="e(n)" symbosz=5 symbol="o" min1=0 max1=%d d1=1'%niter)		
    Graph(init2+'-diffsd','title="Model error 2" label1=Iteration label2="e(n)" symbosz=5 symbol="o" min1=0 max1=%d d1=1'%niter)		

#Making movie
    Plot(init1+'-sigs1',sigs1,'Movie')
    Plot(init2+'-sigs2',sigs2,'Movie')
    Flow(deblended1,sig1,'cp')
    Flow(deblended2,sig2,'cp')

def deblendslet(unblended1,
	  unblended2,
	  blended1,
          blended2,
	  init1,
	  init2,
	  deblended1,
	  deblended2,
	  shottime1,
	  shottime2,
	  mute1,
	  mute2,
	  n1,
	  n2,
	  r1,
	  r2,
	  padno,
	  niter,
	  ddip,
	  mode,
	  thr,
	  pc,
	  fhi,
	  fraction):
    Flow('mask',blended1,'math output=1 | pad n2=%d'%(padno))
    sig1, sig2 = (init1,init2)
    sigs1=[]
    sigs2=[]
    diffsa=[]
    diffsb=[]
    diffsc=[]
    diffsd=[]
    for i in range(niter):
#######################################################################################
# update dip map every "ddip" iterations
#######################################################################################
    	if i % ddip ==0 :
		dip1='dip%d'%int(i/ddip)
		dip2='udip%d'%int(i/ddip)
		Flow(dip1,[sig1,'mask'],
     			'''
			bandpass fhi=%d | pad n2=%d | 
			dip mask=${SOURCES[1]} rect1=%d rect2=%d
			'''%(fhi,padno,r1,r2))
		Flow(dip2,[sig2,'mask'],
     			'''
			bandpass fhi=%d | pad n2=%d | 
			dip mask=${SOURCES[1]} rect1=%d rect2=%d
			'''%(fhi,padno,r1,r2))
#######################################################################################
 	old1,old2 = (sig1,sig2)
	diffa=init1+'-diffa%d'%i
	diffb=init2+'-diffb%d'%i
	diffc=init1+'-diffc%d'%i
	diffd=init2+'-diffd%d'%i

	(nsig1,nsig2)=step(sig1,sig2,blended1,blended2,init1,init2,shottime1,shottime2,i,fraction)
	(sig1,sig2)=shapeslet(nsig1,nsig2,dip1,dip2,mute1,mute2,padno,n2,mode,thr,i)

	Flow(diffa,[old1,sig1],'diff match=${SOURCES[1]} | math output="input/%d/%d"'%(n1,n2))
	Flow(diffb,[old2,sig2],'diff match=${SOURCES[1]} | math output="input/%d/%d"'%(n1,n2))
        Flow(diffc,[unblended1,sig1],'diff match=${SOURCES[1]} | math output="input/%d/%d"'%(n1,n2))
	Flow(diffd,[unblended2,sig2],'diff match=${SOURCES[1]} | math output="input/%d/%d"'%(n1,n2))

    	Greyplot(sig1,'title="Esti R 1 (iter=%d)" pclip=%d'% (i+1,pc))
    	Greyplot(sig2,'title="Esti R 2 (iter=%d)" pclip=%d' % (i+1,pc))	

	sigs1.append(sig1)
    	sigs2.append(sig2)
        diffsa.append(diffa)
        diffsb.append(diffb)
	diffsc.append(diffc)
	diffsd.append(diffd)

    Flow(init1+'-diffsa',diffsa,'cat axis=2 ${SOURCES[1:%d]}|math output="log(input)/log(10)" | transp | put n1=%d o1=1 d1=1 '%(len(diffsa),len(diffsa)))
    Flow(init2+'-diffsb',diffsb,'cat axis=2 ${SOURCES[1:%d]}|math output="log(input)/log(10)" | transp |  put n1=%d o1=1 d1=1'%(len(diffsb),len(diffsb)))
    Flow(init1+'-diffsc',diffsc,'cat axis=2 ${SOURCES[1:%d]}|math output="log(input)/log(10)" | transp |  put n1=%d o1=1 d1=1'%(len(diffsc),len(diffsc)))
    Flow(init2+'-diffsd',diffsd,'cat axis=2 ${SOURCES[1:%d]}|math output="log(input)/log(10)" | transp |  put n1=%d o1=1 d1=1'%(len(diffsd),len(diffsd)))

    Graph(init1+'-diffsa','title="Data error 1" label1=Iteration label2="c(n)" symbosz=5 symbol="o" min1=0 max1=%d d1=1'%niter)	
    Graph(init2+'-diffsb','title="Data error 2" label1=Iteration label2="c(n)" symbosz=5 symbol="o" min1=0 max1=%d d1=1'%niter)
    Graph(init1+'-diffsc','title="Model error 1" label1=Iteration label2="e(n)" symbosz=5 symbol="o" min1=0 max1=%d d1=1'%niter)	
    Graph(init2+'-diffsd','title="Model error 2" label1=Iteration label2="e(n)" symbosz=5 symbol="o" min1=0 max1=%d d1=1'%niter)	

#Making movie
    Plot(init1+'-sigs1',sigs1,'Movie')
    Plot(init2+'-sigs2',sigs2,'Movie')

    Flow(deblended1,sig1,'cp')
    Flow(deblended2,sig2,'cp')
   




