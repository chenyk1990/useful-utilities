#required python package: bruges
#https://github.com/seg/tutorials/blob/master/1606_Wavelet_estimation/Wavelet_estimation_for_well_ties.ipynb
import numpy as np
from matplotlib import pyplot as plt
from bruges.filters import ricker, rotate_phase
from scipy import linalg as la
from numpy.linalg import lstsq

from sklearn import linear_model
#%matplotlib inline

import warnings
# hide all the warnings
warnings.filterwarnings('ignore')

#Read in the data
r = np.load("rpp.npy")
s = np.load("seismic.npy")
dt = .002  # sample rate in seconds
t = np.arange(s.size) * dt

#
fig01 = plt.figure(figsize=(3,8))

# Reflectivity track
ax = fig01.add_subplot(121)
ax.plot(r, t, 'k')
ax.set_xticks([])
ax.set_xlim(-0.5,0.5)
ax.set_ylim(1.5,0)
ax.set_ylabel('two-way time (s)')

# Seismic track
ax2 = fig01.add_subplot(122)
ax2.plot(s, t, 'k')
ax2.fill_betweenx(t, s, 0, s > 0, color='k', alpha=1.0)
ax2.set_ylim(1.5,0)
ax2.set_xlim(-0.25,0.25)
ax2.set_xticks([])
ax2.set_yticks([])
plt.show()

freqs = [5, 80, 130, 160]
c = 1.0
points = c*np.array([-50,-5,-5,-50])



amp_spec = np.abs(np.fft.rfft(s))
f = np.fft.rfftfreq(len(s), d=dt)
P = 20 * np.log10(amp_spec)   #Power in Decibel Scale

fig02 = plt.figure()
ax = fig02.add_subplot(111)
ax.plot(f, P,'m')
ax.plot(freqs, points, 'ko-', lw=2, zorder=2, ms=6)
for fr in freqs:
    ax.text(fr,1.20*points[0], fr, ha='center', va='top', fontsize=15)
ax.set_xlabel('frequency (Hz)', fontsize=14)
ax.set_ylabel('power (dB)', fontsize=14)
plt.show()
# Uncomment this next line if you want to save the figure
# fig02.savefig('figure_1.png', dpi=500)


phase_spec = np.angle(np.fft.rfft(s))
plt.plot(f, phase_spec, 'm')
plt.yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
           [r'$-\pi$', r'$-\pi/2$', r'$0$', r'$+\pi/2$', r'$+\pi$'])
plt.xlabel('frequency (Hz)', fontsize=14)
plt.ylabel('phase', fontsize=14)
plt.show()


from bruges.filters.wavelets import ormsby, ricker
duration = 0.256 # seconds
orms = ormsby(duration, dt=0.002, f=freqs)
tw = np.arange(-len(orms)//2,len(orms)//2,1)*dt

fig03 =plt.figure(figsize=(8,3))
ax = fig03.add_subplot(111)
ax.plot(tw, orms, 'blue', lw=2, alpha=0.75, 
        label='Ormsby (%i-%i-%i-%i Hz)' %(freqs[0],freqs[1],freqs[2],freqs[3]))
ax.legend(loc=1)
ax.set_xlim(-0.064,0.064)
ax.set_xlabel('time (s)', fontsize=14)
ax.set_ylabel('amplitude', fontsize=14)
ax.grid()
plt.show()

# Wavelet estimation by autocorrelation
dw = 64 # number of samples to display on either side of zero
acorr = np.correlate(s, s, mode='same')
w1 = acorr[len(s)//2-dw//2:len(s)//2+dw//2]

def norm(data):
    return data/np.amax(data)
    
fig04 =plt.figure(figsize=(8,3))
ax = fig04.add_axes([0.1, 0.15, 0.8, 0.7])
ax.plot(tw, orms, 'blue', lw=2, alpha=0.75, 
         label='Ormsby (%i-%i-%i-%i Hz)' %(freqs[0],freqs[1],freqs[2],freqs[3]))
ax.plot(np.arange(0,len(w1))*dt-0.064, norm(w1), 'k', lw=2, alpha=0.75,
         label='Autocorrelation')
ax.legend(loc=1)
ax.set_xlim(-0.064,0.064)
ax.set_xlabel('time (s)', fontsize=14)
ax.set_ylabel('amplitude', fontsize=14)
ax.grid()
fig04.tight_layout

# Uncomment this next line if you want to save the figure
# fig04.savefig('figure_2.png', dpi=500)


#Wavelet estimation by spectral division

def spectral_division(reflectivity, data):
    
    seis_fft = np.fft.fft(data)
    ref_fft = np.fft.fft(reflectivity)

    wavelet_spec = seis_fft / ref_fft
    wavelet_div = np.fft.ifft(wavelet_spec)
    
    return wavelet_div
        

spec_div = spectral_division(r, s)

fig05 = plt.figure()
ax = fig05.add_subplot(111)
ax.plot(t, spec_div, 'k', lw=2)
ax.set_xlim([0,0.256])
plt.show()

def wigner(rpp, seismic):
    opConvolve = la.toeplitz(rpp)
    wavelet = lstsq(opConvolve, seismic)[0]
    return wavelet
    
wigner_wave = wigner(r, s)


fig06 = plt.figure()
ax = fig06.add_subplot(111)
ax.plot(t, wigner_wave, 'k', lw=2)
ax.set_xlim([0,0.128])
ax.set_ylim([-1,1])    


Y1 = 20 * np.log10(np.abs(np.fft.rfft(s))) 
R1 = 20 * np.log10(np.abs(np.fft.rfft(r))) 
W1 = 20 * np.log10(np.abs(np.abs(np.fft.rfft(s) / np.fft.rfft(r)))) 
W2 = 20 * np.log10(np.fft.rfft(wigner_wave)) 


fig07 = plt.figure(figsize=(15,3))
ax = fig07.add_axes([0.1, 0.15, 0.25, 0.7])
ax.plot(f,Y1, 'k', lw=1)
ax.set_title('data', fontsize=14)
ax.set_ylabel('power (dB)', fontsize=14)
ax.set_xlabel('frequency (Hz)', fontsize=14)
ax.set_ylim(-80,10)

ax2 = fig07.add_axes([0.1 + 1*0.8/3, 0.15, 0.25, 0.7])
ax2.plot(f,R1, 'k', lw=1)
ax2.set_title('reflectivity', fontsize=14)
ax2.set_xlabel('frequency (Hz)', fontsize=14)
ax2.set_ylim(-80,10)
ax2.set_yticklabels([])

ax3 = fig07.add_axes([0.1 + 2*0.8/3, 0.15, 0.25, 0.7])
ax3.plot(f,W1, 'k', lw=1)
# ax3.plot(X,W2, 'dark blue', lw=1)
ax3.set_title(' "wavelet" ', fontsize=14)
ax3.set_xlabel('frequency (Hz)', fontsize=14)
ax3.set_ylim(-80,10)
ax3.set_yticklabels([])
plt.show()

# Uncomment this next line if you want to save the figure
# fig07.savefig('figure_3.png', dpi=500)

#Wavelet estimation by least squares

clf = linear_model.Ridge(alpha = 0.5, fit_intercept=False)
R = la.toeplitz(r)
clf.fit(R, s)
wavelet = clf.coef_

Y2 = 20* np.log10(np.abs(np.fft.rfft(wavelet)))  

fig08 = plt.figure(figsize = (10,3))
ax = fig08.add_subplot(121)
ax.plot(t, wavelet, 'k', lw=2)
ax.set_xlim([0,0.128])
ax.set_title('wavelet', fontsize=14)
ax.set_ylabel('amplitude', fontsize=14)
ax.set_xlabel('time (s)', fontsize=14)

# Check the spectra
ax2 = fig08.add_subplot(122)
ax2.plot(f,Y2, 'k', lw=2)
ax2.set_title('spectrum of wavelet', fontsize=14)
ax2.set_ylabel('power (dB)', fontsize=14)
ax2.set_xlabel('frequency (Hz)', fontsize=14)
plt.show()

# modelled seismic
fig09 = plt.figure(figsize=(15,3))
ax = fig09.add_subplot(111) 
ax.plot(t, np.dot(R, wavelet), 'k', lw=2)
ax.set_title('synthetic', fontsize=14)
ax.set_ylabel('amplitude', fontsize=14)
ax.set_xlabel('time (s)', fontsize=14)


wavelet_size = 15 # samples
opProj = np.zeros((r.size, r.size))
opProj[:wavelet_size, :wavelet_size] = np.eye(wavelet_size)

op  = np.dot(R, opProj)
wavelet = lstsq(op, s)[0]

Y3 = 20* np.log10(np.abs(np.fft.rfft(wavelet))) 

fig10 = plt.figure(figsize = (10,3))
ax = fig10.add_axes([0.1, 0.15, 0.3, 0.7])
ax.plot(t, wavelet, 'k', lw=2)
ax.set_xlim([0,0.128])
ax.set_title('wavelet', fontsize=14)
ax.set_ylabel('amplitude', fontsize=14)
ax.set_xlabel('time (s)', fontsize=14)

# Check the spectra
ax2 = fig10.add_axes([0.5, 0.15, 0.35, 0.7])
ax2.plot(f, Y3, 'k', lw=2)
ax2.set_title('spectrum of wavelet', fontsize=14)
ax2.set_ylabel('power (dB)', fontsize=14)
ax2.set_xlabel('frequency (Hz)', fontsize=14)

# Uncomment this next line if you want to save the figure
# fig10.savefig('figure_4.png', dpi=500)


#Create synthetics
synth = np.dot(op, wavelet)

# Params
ylim = (0.75, 1.5)
# N wiggles
nwigs = 5
gap = 0.1
g = 0.25  # gain

fig10 = plt.figure(figsize=(6,8))
ax1 = fig10.add_subplot(141)
ax2 = fig10.add_subplot(142)
ax3 = fig10.add_subplot(143)

# Reflectivity track
ax1.plot(r, t, 'b')
ax1.set_xticks([])
ax1.set_xlim(-0.5,0.5)
ax1.set_xticks([])
ax1.set_xlabel('reflectivity', fontsize=14)
ax1.xaxis.set_label_position('top') 
ax1.set_ylim(ylim[1],ylim[0])
ax1.set_ylabel('two-way time (s)', fontsize=14)

# Synthetic track
for i in range(nwigs):
    ax2.plot(synth+i*gap, t, 'k', lw=0.5)
    ax2.fill_betweenx(t, synth+i*gap, i*gap,
                      synth+i*gap > i*gap, color='k', alpha=1.0)
ax2.set_xlim(-gap,gap*(nwigs))
ax2.set_ylim(ylim[1],ylim[0])
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_xlabel('synthetic', fontsize=14)
ax2.xaxis.set_label_position('top') 


# Seismic track
for i in range(nwigs):
    ax3.plot(s+i*gap, t, 'k', lw=0.5)
    ax3.fill_betweenx(t, s+i*gap, i*gap,
                      s+i*gap > i*gap, color='k', alpha=1.0)
ax3.set_xlim(-gap,gap*(nwigs))
ax3.set_ylim(ylim[1],ylim[0])
ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_xlabel('data', fontsize=14)
ax3.xaxis.set_label_position('top') 

# Wavelet track
if False:
    ax4 = fig.add_subplot(144)
    ax4.plot(wavelet, t, 'k', lw=0.5)
    ax4.fill_betweenx(t, wavelet, 0,
                      wavelet > 0, color='k', alpha=1.0)
    ax4.set_ylim(ylim[1],ylim[0])
    ax4.set_xlim(-0.1,0.1)
    ax4.set_xticks([])
    ax4.set_yticks([])
    ax4.tight_layout()

# Uncomment this next line if you want to save the figure
# fig10.savefig('figure_5.png', dpi=500)
plt.show()

np.corrcoef(s,synth)[0][1]

wzeroph = norm(np.correlate(norm(wavelet), norm(wavelet), mode='same'))
dw = 64
wshort = wzeroph[len(wzeroph)//2-dw/2:len(wzeroph)//2+dw/2]
tshort = np.arange(0,len(wshort),1)*dt

plt.plot(tshort, wshort,'k', lw=2)
plt.fill_between(tshort, wshort, -0.005,
                      wshort > -0.005, color='k', alpha=0.85)
plt.show()


fig11 =plt.figure(figsize=(8,3))
ax = fig11.add_subplot(111)
ax.plot(tshort-0.064, wshort,'k', lw=2, alpha=0.5, label='estimated wavelet')

ax.plot(tw, orms, 'blue', lw=2, alpha=0.5, 
         label='Ormsby (%i-%i-%i-%i Hz)' %(freqs[0],freqs[1],freqs[2],freqs[3]))
ax.legend(loc=1)
ax.set_xlim(-0.064,0.064)
ax.set_xlabel('time (s)', fontsize=14)
ax.set_ylabel('amplitude', fontsize=14)
ax.grid()
plt.show()

synth = np.convolve(norm(wzeroph),r,mode='same')
# Params
ylim = (0.75, 1.5)
# N wiggles
nwigs = 5
gap = 0.1
g = 0.25  # gain

fig12 = plt.figure(figsize=(6,8))
ax1 = fig12.add_subplot(141)
ax2 = fig12.add_subplot(142)
ax3 = fig12.add_subplot(143)

# Reflectivity track
ax1.plot(r, t, 'b')
ax1.set_xticks([])
ax1.set_xlim(-0.5,0.5)
ax1.set_xticks([])
ax1.set_xlabel('reflectivity', fontsize=14)
ax1.xaxis.set_label_position('top') 
ax1.set_ylim(ylim[1],ylim[0])
ax1.set_ylabel('two-way time (s)', fontsize=14)

# Synthetic track
for i in range(nwigs):
    ax2.plot(g*synth+i*gap, t, 'k', lw=0.5)
    ax2.fill_betweenx(t, g*synth+i*gap, i*gap,
                      g*synth+i*gap > i*gap, color='k', alpha=1.0)
ax2.set_xlim(-gap,gap*(nwigs))
ax2.set_ylim(ylim[1],ylim[0])
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_xlabel('synthetic', fontsize=14)
ax2.xaxis.set_label_position('top') 


# Seismic track
for i in range(nwigs):
    ax3.plot(s+i*gap, t, 'k', lw=0.5)
    ax3.fill_betweenx(t, s+i*gap, i*gap,
                      s+i*gap > i*gap, color='k', alpha=1.0)
ax3.set_xlim(-gap,gap*(nwigs))
ax3.set_ylim(ylim[1],ylim[0])
ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_xlabel('data', fontsize=14)
ax3.xaxis.set_label_position('top') 

# Wavelet track
if False:
    ax4 = fig12.add_subplot(144)
    ax4.plot(wavelet, t, 'k', lw=0.5)
    ax4.fill_betweenx(t, wavelet, 0,
                      wavelet > 0, color='k', alpha=1.0)
    ax4.set_ylim(ylim[1],ylim[0])
    ax4.set_xlim(-0.1,0.1)
    ax4.set_xticks([])
    ax4.set_yticks([])
    ax4.tight_layout()
plt.show()
np.corrcoef(s,synth)[0][1]
    
                      
                      



    

















































