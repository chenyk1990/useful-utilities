#ipynb file path:
#specf/seg/1604_Function_of_interpolation

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d



# Make a 'trace' of amplitudes.
amps = np.array([0,1,3,0,-2,-5,-2,2,5,0])

# Give the trace a time basis.
time = np.linspace(0, 18, amps.size)
print(time)


f = interp1d(time, amps, kind='linear')
print(f(3.1415926))


x, y = 3.14159, f(3.14159)
plt.figure(figsize=(15, 3))
plt.plot(time, amps, 'o-')
plt.axhline(0, c='k')

# Plot our point.
plt.plot(x, y, 'ro')
plt.axvline(x, ymin=0, ymax=(6+y)/12, c='r')
plt.axhline(y, xmin=0, xmax=x/18, c='r')
plt.show()

#Different measures of amplitude at a given
trace = np.array([0,1,3,0,-2,-5,-2,2,5,0,-1,1,0,3,3,-1,-3,-5,-1,0])
horizon = 7.544  # No units, it's just a sort of fractional index.

def read_amp(trace, time, method='linear'):
    """Wraps ``interp1d``.
    """
    indices = np.arange(trace.size)
    
    # Create the interpolation function, f.
    f = interp1d(indices, trace, kind=method, assume_sorted=True)
    
    return f(time)
    
hor_spline = read_amp(trace, horizon, method='cubic')
print(hor_spline)


hor_linear = read_amp(trace, horizon, method='linear')
print(hor_linear)


zero_order = read_amp(trace, horizon, method='zero')
print(zero_order)


# Set up basis vectors
x = np.arange(trace.size)
t = np.linspace(0, trace.size-1, 2000)
amp_spline = read_amp(trace, t, method='cubic')
amp_linear = read_amp(trace, t, method='linear')
amp_nearest = read_amp(trace, t, method='nearest')

styles = {'ha':'left', 'va':'center', 'size':14}
fig = plt.figure(figsize=(16,4))

#ax = fig.add_axes([0.05,0.15,0.825,0.75])
ax = fig.add_subplot(111)

# Fix some types
hor_spline = float(hor_spline)
hor_linear = float(hor_linear)
zero_order = float(zero_order)

# Samples
ax.scatter(x, trace, s=50, c='k', lw=0, zorder=10)
for i, v in enumerate(trace):
    adj = 0.01
    if v > 0:
        ax.axvline(x=i, ymin=6/13, ymax=6/13 + v/13., c='k', alpha=0.67, lw=2, zorder=9)
    else:
        ax.axvline(x=i, ymin=6/13 + v/13., ymax=6/13-adj, c='k', alpha=0.67, lw=2, zorder=9)

# Spline interpolated
ax.plot(t, amp_spline, 'k', lw=3)
ax.fill_between(t, 0, amp_spline, where=amp_spline>0, color='k', alpha=0.5, lw=0)
ax.text(0.3, 2.7, "spline", color='k', **styles)

# Nearest interpolated
ax.plot(t, amp_nearest, '#22bb22', lw=2)
ax.fill_between(t, 0, amp_nearest, where=amp_nearest>0, color='#22bb22', alpha=0.5, lw=0)
ax.text(2.3, 3.7, "nearest", color='#22bb22', **styles)

# Linear interpolated
ax.plot(t, amp_linear, '#00ccff', lw=2)
ax.fill_between(t, 0, amp_linear, where=amp_linear>0, color='#00aaff', alpha=0.5, lw=0)
ax.text(2.8, 1.6, "linear", color='#00aaff', **styles)

# Horizon crossing point and amplitudes
pos = 10.1
ax.axvline(x=horizon, c='r', lw=2)
ax.text(horizon+0.2, -5, "HORIZON at t = {}".format(horizon), color='r', **styles)

styles = {'ha':'left', 'va':'center', 'size':12}
ax.axhline(y=5, xmin=8/19., xmax=10/19, c='r', lw=1)
ax.axhline(y=hor_spline, xmin=horizon/19., xmax=10/19, c='r', lw=1)
ax.axhline(y=hor_linear, xmin=horizon/19., xmax=10/19, c='r', lw=1)
#ax.axhline(y=zero_order, xmin=7/19., xmax=10/19, c='r', lw=1)
ax.text(pos, 5+0.4, "nearest = {:.2f}".format(5), color='#22bb22', **styles)
ax.text(pos, hor_spline, "spline = {:.2f}".format(hor_spline), color='k', **styles)
ax.text(pos, hor_linear-0.3, "linear = {:.2f}".format(hor_linear), color='#00aaff', **styles)
#ax.text(pos, zero_order, "zero-order = {:.2f}".format(zero_order), color='k', **styles)

# Axes etc.
ax.text(0.3, 5.5, 'Sample interpolation', size=15, weight="bold")
ax.axhline(0, color='k')
ax.set_xlim(0, 19)
ax.set_ylim(-6, 7)
plt.xlabel('time samples', size=13, ha='right')
ax.xaxis.set_label_coords(1.0, -0.025)
plt.ylabel('amplitude', size=13)
plt.minorticks_on()
ax.tick_params(length=4, width=1, which='both')
ax.set_yscale('linear', subsx=np.arange(-6, 6, 2))  
plt.gca().xaxis.grid(True, which='both', lw=2, color='k', alpha=0.1)

# Save figure to your home directory.
plt.savefig("trace_sampling.png", dpi=400)

# Show figure.
plt.show()

#Real data
#data path
#https://s3.amazonaws.com/agilegeo/Penobscot.h5 

import h5py
h5f = h5py.File('data/Penobscot.h5','r')
seismic = h5f['amplitude'][:]
h5f.close()


clip = np.percentile(seismic, 99)
fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(111)
plt.imshow(seismic[:,100,:], cmap="Greys", vmin=-clip, vmax=clip)
plt.colorbar(label="Amplitude", shrink=0.8)
ax.set_xlabel("Trace number")
ax.set_ylabel("Time sample")
plt.show()

#Load the horizon
horizon = np.load('data/Penobscot_Seabed.npy')

def regularize_horizon(horizon, extents, adj=(0,0)):
    output = np.empty(extents)
    output[:] = np.nan
    adj_x, adj_y = adj
    for (x, y, z) in horizon:
        output[int(x+adj_x), int(y+adj_y)] = z
    return output
    

#Check the horizon looks OK
fig = plt.figure(figsize=(15,6))
ax = fig.add_subplot(111)
plt.imshow(-horizon, aspect=0.5, cmap="gist_earth", origin='lower')
plt.colorbar(label="Two-way time")
ax.set_xlabel("Trace number")
ax.set_ylabel("Trace number")
plt.show()    


#Extract amplitude on one line
line = 300

horizon[horizon==0] = np.nan

from matplotlib.font_manager import FontProperties

fig = plt.figure(figsize=(18,8))
ax = fig.add_subplot(111)
plt.imshow(seismic[:,line,:], cmap="Greys", vmin=-clip, vmax=clip)
plt.plot(horizon[line], label='Seabed', c='r', lw=2)

# Colourbar
cbar = plt.colorbar(label="Amplitude", shrink=0.9)
cbar.ax.tick_params(labelsize=12) 
text = cbar.ax.yaxis.label
font = FontProperties(size=14)
text.set_font_properties(font)

# Axes etc.
plt.legend(fontsize=15)
plt.xlim(0, 480)
plt.ylim(250, 0)
ax.set_xlabel("Trace number", size=14)
ax.set_ylabel("Time sample", size=14)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.show()

fig = plt.figure(figsize=(18,6))
ax = fig.add_subplot(111)
plt.imshow(seismic[25:76,line,:151], cmap="Greys", vmin=-clip, vmax=clip, interpolation='none')
ax.plot(horizon[line]-25, label='Seabed', c='r', lw=2)

# Vertical line
pos = 50
twt = horizon[line][pos]
ax.axvline(x=pos, c='#00ddff', lw=2)
ax.text(pos+1, 4, "trace number {}".format(pos), color='#00eeff', fontsize=14, fontweight='bold')
ax.text(pos+1, 7, "TWTT {:.2f} samples".format(twt), color='#ee1111', fontsize=14, fontweight='bold')

# Colourbar
# cbar = plt.colorbar(label="Amplitude", shrink=0.8)
# cbar.ax.tick_params(labelsize=12) 
# text = cbar.ax.yaxis.label
# font = FontProperties(size=14)
# text.set_font_properties(font)


# Axes etc.
plt.legend(fontsize=15)
plt.xlim(0, 150)
plt.ylim(50, 0)
ax.set_xlabel("Trace number", size=14)
ax.set_ylabel("Time sample", size=14)
plt.xticks(fontsize = 12) 
plt.yticks([0, 10, 20, 30, 40, 50], [25, 35, 45, 55, 65, 75])
plt.yticks(fontsize = 12) 

plt.savefig("seismic_horizon_segment.png", dpi=250)
plt.show()

def get_amp(seismic, horizon, method="nearest"):
    # Deal with NaNs in the input and prepare the output.
    horizon[np.isnan(horizon)] = 0
    amps = np.zeros_like(horizon)
    amps[:] = np.nan
    
    # Iterate over traces and get the interpolated value.
    if seismic.ndim == 1:
        return read_amp(seismic, horizon, method)
    elif seismic.ndim == 2:
        # Treat like 2D
        for j, (trace, t) in enumerate(zip(seismic.T, horizon)):
            amps[j] = read_amp(trace, t, method)
    elif seismic.ndim == 3:
        # Treat like 3D
        for i, section in enumerate(np.swapaxes(seismic.T, 0, 1)):
            for j, (trace, t) in enumerate(zip(section, horizon[i])):
                amps[i, j] = read_amp(trace, t, method)
    else:
        raise Exception("Too many dimensions.")
 
    amps[amps==0] = np.nan
    return amps


a = get_amp(seismic[:, line, :], horizon[line], method="nearest")
fig = plt.figure(figsize=(18,4))
ax = fig.add_subplot(111)
ax.plot(a)
ax.set_xlabel("Trace number")
ax.set_ylabel("Amplitude")
plt.show()

#Extract amplitude on entire horizon

amp_full = get_amp(seismic, horizon, method="linear")
    

fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(111)

# Draw
plt.imshow(amp_full, cmap='Greys', origin='lower', aspect=0.5, vmin=0)
ax.set_ylabel("Trace number", size=14)
ax.set_xlabel("Trace number", size=14)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 

# Colourbar
cbar = plt.colorbar(label="Amplitude", shrink=0.9)
cbar.ax.tick_params(labelsize=12) 
text = cbar.ax.yaxis.label
font = FontProperties(size=14)
text.set_font_properties(font)

# Save
plt.savefig("horizon_slice.png", dpi=250, transparent=True)
plt.show()


methods = ["cubic", "quadratic", "linear", "nearest", "zero"]

results = {}
for method in methods:
    results[method] = get_amp(seismic[:,line,:], horizon[line], method=method)

s1, s2 = 0, 151
fig = plt.figure(figsize=(18, 5))
ax = fig.add_subplot(111)
for method, data in results.items():
    ax.plot(data[s1:s2], label=method)
ax.set_xlabel("Trace number in this data segment")
ax.set_ylabel("Amplitude")
ax.legend(loc='lower right')
plt.grid()
plt.show()


#Figure for print

def rms(a):
    mean_squares = np.nansum(a**2.)/a.size
    return np.sqrt(mean_squares)

# Deal with edge
for m, d in results.items():
    d[d<1] = np.nan

params = {'linear': {'c': 'darkturquoise',
                     'lw': 1.5,
                     'zorder': 1,
                     'label': 'linear, 40 ms, RMS {:.0f}'
                    },
          'nearest': {'c': 'limegreen',
                     'lw': 1.5,
                     'zorder': 1,
                     'label': 'nearest, 33.7 ms, RMS {:.0f}'
                    },
          'cubic': {'c': 'k',
                     'lw': 3.0,
                     'zorder': 0,
                     'label': 'cubic, 10 900 ms, RMS 0'
                    },
         }

s1, s2 = 0, 151
r = results.copy()
c = r.pop('cubic')

fig = plt.figure(figsize=(18, 4))
ax = fig.add_subplot(111)
for m, data in results.items():
    if m in ['zero', 'quadratic']: continue
    diff = c - data
    ax.plot(data[s1:s2], label=params[m].pop('label').format(rms(diff)), **params[m])
ax.legend(loc='lower right', fontsize=15)
ax.set_ylabel("Amplitude", size=14)
ax.xaxis.tick_top()
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12)
plt.xlim(0, 150)
plt.grid()

plt.savefig("amplitude_extraction.png", dpi=250)
plt.show()


#What's the difference?
#for method in methods:
#    print(method, end=',')
#    %timeit get_amp(seismic[:,line,:], horizon[line], method=method)
    
def rms(a):
    mean_squares = np.nansum(a**2.)/a.size
    return np.sqrt(mean_squares)
    
    
print("Mean RMS difference between cubic result and...")
r = results.copy()
c = r.pop('cubic')
for method, data in r.items():
    diff = c - data
    print("{}: {:.2f}".format(method, rms(diff)))
    

timeslice = np.zeros_like(horizon) + 124.5


results = {}
for method in methods:
    results[method] = get_amp(seismic[:,line,:], timeslice[line], method=method)
    
    
s1, s2 = 100, 200
fig = plt.figure(figsize=(18, 8))
ax = fig.add_subplot(111)
for method, data in results.items():
    plt.plot(data[s1:s2], label=method)
ax.set_xlabel("Trace number in this data segment")
ax.set_ylabel("Amplitude")
plt.legend()
plt.grid()
plt.show()       


print("Mean RMS difference between cubic result and...")
r = results.copy()
c = r.pop('cubic')
for method, data in r.items():
    diff = c - data
    print("{}: {:.2f}".format(method, rms(diff)))


#Spatial interpolation

from scipy.interpolate import splprep, splev
from mpl_toolkits.mplot3d import Axes3D

trajectory = np.array([[   0,   0,    0],
                       [   0,   0, -100],
                       [   0,   0, -200],
                       [   5,   0, -300],
                       [  10,  10, -400],
                       [  20,  20, -500],
                       [  40,  80, -650],
                       [ 160, 160, -700],
                       [ 600, 400, -800],
                       [1500, 960, -800]])

# Shift the 'well' to its tophole location.
trajectory[:,0] += 3000
trajectory[:,1] += 3000

# Fit the spline and fetch 400 points.
knees, _ = splprep(trajectory.T, s=3.0, k=3)
spline = splev(np.linspace(0, 1, 1000), knees)         
plt.figure(figsize=(12,7))
plt.gca(projection='3d')
plt.plot(*spline, color='grey', lw=3, alpha=0.75)
plt.show()

path = np.array(spline).T[::-1]

seismic.shape

s = np.swapaxes(seismic, 0, 2)
s.shape


from scipy.interpolate import RegularGridInterpolator

# Make linear spaces in real-world units for the seismic domain.
x = np.linspace(0, 12000, 481) # m
y = np.linspace(0, 7500, 601)  # m
t = np.linspace(-1000, 0, 251)  # ms

f = RegularGridInterpolator((x, y, t), s)

amp_along_path = f(path)


fig = plt.figure(figsize=(18,3))
ax = fig.add_subplot(111)
x = np.arange(1000)
ax.plot(amp_along_path, 'k')
ax.fill_between(x, 0, amp_along_path, where=amp_along_path>0, color='k', alpha=0.75, lw=0)
plt.show()

f([123, 456, -789])

# We need the data in yet another order!
s = np.swapaxes(s, 0, 1)

x = np.linspace(1000, 1600, 601)  # Inline numbers
y = np.linspace(1000, 1480, 481)  # Crossline numbers
t = np.linspace(0, 1000, 251)  # ms

g = RegularGridInterpolator((x, y, t), s, method="nearest")  # or 'linear' (default)

horizon_data = np.loadtxt('data/Penobscot_Seabed.txt')
horizon_data

amplitudes = np.copy(horizon_data)
amplitudes[:,2] = g(horizon_data)

amps = regularize_horizon(amplitudes, (601, 481), adj=(-1000, -1000))


fig = plt.figure(figsize=(16,8))
ax = fig.add_subplot(111)
plt.imshow(amps, cmap='Greys', origin='lower', aspect=0.5, vmin=0)
ax.set_ylabel("Trace number")
ax.set_xlabel("Trace number")
plt.colorbar(label="Amplitude")
plt.show()










































