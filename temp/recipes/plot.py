## This is a script for the ploting different type of figures using a template.

## Copyright (C) 2013
## Yangkang Chen, The Univeristy of Texas at Austin


from rsf.proj import*
def Grey(data,other): 
	Result(data,'grey label2=Trace unit2="" label1=Time unit1="s" title="" labelsz=10 labelfat=4 font=2 titlesz=10 titlefat=4 screenht=10.24 screenratio=0.75 screenwd=13.65 wherexlabel=b wheretitle=t color=b bartype=v %s'%other)

def Greyplot(data,other): 
	Plot(data,'grey label2=Trace unit2="" label1=Time unit1="s" title="" labelsz=10 labelfat=4 font=2 screenht=10.24 screenratio=0.75 screenwd=13.65 wherexlabel=b wheretitle=t color=b bartype=v %s'%other)

def Graph(data,other):
	Result(data,'graph label1="" label2="" unit1="" unit2=""  title="" labelsz=10 labelfat=4 font=2 screenht=10.24 screenratio=0.75 screenwd=13.65 wherexlabel=b wheretitle=t %s' %other)

def Graphplot(data,other):
	Plot(data,'graph label1="" label2="" unit1="" unit2=""  title="" labelsz=10 labelfat=4 font=2 screenht=10.24 screenratio=0.75 screenwd=13.65 wherexlabel=b wheretitle=t %s' %other)

def Wig(data,other):
	Result(data,'window j2=4 | wiggle transp=y yreverse=y poly=y clip=0.1 labelsz=10 labelfat=4 font=2 screenht=10.24 screenratio=0.75 screenwd=13.65 wanttitle=n')
   
def Wigplot(data,other):
	Plot(data,'window j2=4 | wiggle transp=y yreverse=y poly=y clip=0.1 labelsz=10 labelfat=4 font=2 screenht=10.24 screenratio=0.75 screenwd=13.65 wanttitle=n')

def Dots(data,other):
	Result(data,'dots %s'%other)

def Dotsplot(data,other):
	Plot(data,'dots %s'%other)

def Grey3(data,other):
	Result(data,
       '''
       byte clip=0.0001  |
       transp plane=23 |
       grey3 flat=n  frame1=500 frame2=125 frame3=5
       title=Data point1=0.8 point2=0.8  %s pclip=5
       '''%other)
#

def Grey3(data,other):
	Result(data,
		'''
		byte | grey3 flat=n label1=Time wanttitle=n labelsz=10
		unit1=s label3=Inline label2=Crossline
		point1=0.8 point2=0.8 frame1=80 frame2=60 frame3=80 %s
		'''%other)

def Vel(data,other):
	Result(data,
     '''
     grey color=j allpos=y bias=1.5 clip=0.7
     scalebar=y barreverse=y barunit=km/s
     label2=Midpoint unit2=km label1=Time unit1=s
     title="NMO Velocity"  %s
     '''%other )

def Pick(data,other):
	 Result(data,
       '''
       byte allpos=y gainpanel=all |
       transp plane=23 |
       grey3 flat=n frame1=500 frame2=125 frame3=25 
       label1=Time unit1=s color=j
       label3=Velocity unit3=km/s 
       label2=Midpoint unit2=km
       title="Velocity Scan" point1=0.8 point2=0.8 %s
       '''%other)

def Pickmovie(data,other):
	 Result(data,
       '''
       byte allpos=y gainpanel=all |
       transp plane=23 |
       grey3 flat=n frame1=500 frame2=125 frame3=25 
       label1=Time unit1=s color=j
       label3=Velocity unit3=km/s 
       label2=Midpoint unit2=km movie=2 dframe=5
       title="Velocity Scan" point1=0.8 point2=0.8 %s
       '''%other)

## Create label A
Plot('labela',None,
	'''
	box x0=9 y0=7.55 label="A" xt=0.5 yt=-0.5 length=0.75 
	''')

## Create label B
Plot('labelb',None,
	'''
	box x0=7 y0=4.6 label="B" xt=-0.5 yt=0.5 length=0.75
	''')



