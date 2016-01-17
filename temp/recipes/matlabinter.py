from rsf.proj import*
from rsf.prog import RSFROOT


########################################################################
# INITIALIZATION
########################################################################
matlab         = WhereIs('matlab')
matROOT = '../matfun'
matfun = 'Function'
matlabpath = os.environ.get('MATLABPATH',os.path.join(RSFROOT,'lib'))

if not matlab:
    sys.stderr.write('\nCannot find Matlab.\n')
    sys.exit(1)

############################################################
## no parameter
############################################################
Flow('file1 file2',[os.path.join(matROOT,matfun+'.m'),'file0'],
     '''MATLABPATH=%(matlabpath)s %(matlab)s 
     -nosplash -nojvm -r "addpath %(matROOT)s;%(matfun)s('${SOURCES[1]}','${TARGETS[0]}','${TARGETS[1]}');quit"
     '''%vars(),stdin=0,stdout=-1)


value1=1
value2=0.1

############################################################
## with parameter
############################################################
Flow('file1',[os.path.join(matROOT,matfun+'.m'),'file0'],
     '''MATLABPATH=%(matlabpath)s %(matlab)s 
     -nosplash -nojvm -r "addpath %(matROOT)s;%(matfun)s('${SOURCES[1]}',%(value1)d,%(value2)g,'${TARGETS[0]}');quit"
     '''%vars(),stdin=0,stdout=-1)


###Using FXEMD to enhance useful horizontal signals
########################################################################
# INITIALIZATION
########################################################################
matlab         = WhereIs('matlab')
matROOT = '/home/chenyk/chenyk/matlibcyk/Mada_matlab'
matfun = 'FXEMD'
matlabpath = os.environ.get('MATLABPATH',os.path.join(RSFROOT,'lib'))

if not matlab:
    sys.stderr.write('\nCannot find Matlab.\n')
    sys.exit(1)
n1=
n2=
dt=0.004
lf=5
hf=120
N=1
verb=0
############################################################
## with parameter
############################################################
Flow('flat-emd0',[os.path.join(matROOT,matfun+'.m'),'flat'],
     '''MATLABPATH=%(matlabpath)s %(matlab)s 
     -nosplash -nojvm -r "addpath %(matROOT)s;%(matfun)s('${SOURCES[1]}','${TARGETS[0]}',%(n1)d,%(n2)d,%(dt)g,%(lf)g,%(hf)g,%(N)d,%(verb)d);quit"
     '''%vars(),stdin=0,stdout=-1)
Flow('file-emd','file-emd0','put d2= d1= o2= o1=')

