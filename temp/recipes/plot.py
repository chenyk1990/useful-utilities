## This is a script for the ploting different type of figures using a template.

## Copyright (C) 2013
## Yangkang Chen, The Univeristy of Texas at Austin


from rsf.proj import*
def Grey(data,other): 
	Result(data,'grey label2=Trace unit2="" label1=Time unit1="s" title="" labelsz=10 labelfat=4 font=2 titlesz=10 titlefat=4 screenht=10.24 screenratio=0.75 screenwd=13.65 wherexlabel=b wheretitle=t color=b bartype=v %s'%other) 
	
#font=104

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
       byte allpos=y gainpanel=all bar=bar.rsf |
       transp plane=23 |
       grey3 flat=n frame1=500 frame2=125 frame3=25 
       label1=Time unit1=s color=j
       label3=Velocity unit3=km/s 
       label2=Midpoint unit2=km
       title="Velocity Scan" point1=0.8 point2=0.8 scalebar=y %s
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

## Movie
Result('gxy3-ts',tmp, 'Movie')

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

## Creating framebox
x=9   #spatial
y=1 #vertical
w=2.5  #spatial
w1=1#vertical

Flow('frame.asc',None,'echo %s n1=10 data_format=ascii_float in=$TARGET'% \
	string.join(map(str,(x,y,x+w,y,x+w,y+w1,x,y+w1,x,y))))
Plot('frame','frame.asc',
	'''
	dd type=complex form=native |
	graph min1=8.0 max1=16.0 min2=0 max2=2.5 pad=n plotfat=15 plotcol=2 
	wantaxis=n wanttitle=n yreverse=y 
	''')

## Latex usage	
\bibliographystyle{seg}
\bibliography{bibname}
\newpage
\listoffigures
\newpage
\AtEndDocument{
\begin{figure}
	\centering
	\subfloat[]{\includegraphics[width=0.55\textwidth]{Fig/p1}}
    \subfloat[]{\includegraphics[width=0.55\textwidth]{Fig/p2}}
    \subfloat[]{\includegraphics[width=0.55\textwidth]{Fig/p3}}
	\caption{This is caption.}
	\label{fig:p1,p2,p3}
\end{figure}
}

\begin{table}[h]
\caption{Comparison of computational cost of both the MSSA and the proposed HRSC framework for different data sizes.}
\begin{center}
     \begin{tabular}{|c|c|c|c|} 
	  \hline Test & 300$\ast$20$\ast$20  &  300$\ast$40$\ast$20  &  300$\ast$80$\ast$20\\ 
	  \hline MSSA (sec) & 207.46 & 597.60 & 2473.50 \\
      	  \hline Proposed HRSC (sec)& 	131.78 & 396.89 & 1706.31		       		 \\ 
          \hline
    \end{tabular} 
\end{center}
\label{tbl:compucost}
\end{table}

##basemap
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
plt.figure(figsize=(8,8)) 
# make sure the value of resolution is a lowercase L,
#  for 'low', not a numeral 1
m = Basemap(projection='merc', lat_0=38, lon_0=138,
    resolution = 'h', area_thresh = 1000.0,
    llcrnrlon=125, llcrnrlat=30,
    urcrnrlon=146, urcrnrlat=47)
 
m.drawcoastlines()
m.drawcountries()
m.drawstates()
m.fillcontinents()
m.drawparallels(np.arange(-90., 120.,10.),labels=[1,1,0,0])
m.drawmeridians(np.arange(0., 420., 10.),labels=[0,0,1,1])
m.drawmapboundary()

lon = 135
lat = 37
xm,ym = m(lon, lat)
m.plot(xm, ym, 'bo', markersize=12)
plt.text(xm,ym,'AM',fontsize=20,fontweight='bold',
                    ha='center',va='center',color='r')
plt.show()

#matlab plot
  ylabel('Time(s)','Fontsize',16,'fontweight','bold');
  xlabel('Trace','Fontsize',16,'fontweight','bold');
  title('Noise','Fontsize',16,'fontweight','bold');
  set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold'); 
  print(gcf,'-depsc','-r200','fig.eps');  
  
  

