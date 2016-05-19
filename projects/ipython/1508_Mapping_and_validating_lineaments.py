#https://github.com/seg/tutorials/blob/master/1508_Mapping_and_validating_lineaments/1508_Mapping_and_validating_lineaments.ipynb


import numpy as np
import scipy as sp
from fatiando.gravmag import transform
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as clr
from pylab import imread, imshow, gray, mean
from skimage.filters import threshold_otsu, rank, threshold_adaptive    # N.B.depending on the version of Python you are running
                                                                        # this may have to be "from skimage.filter import ..."
    
from skimage.morphology import disk, closing, opening, erosion, dilation, reconstruction, skeletonize, remove_small_objects
from skimage import img_as_uint, img_as_float
from skimage import io
from skimage import color


# open colormap file, import raw RGB values
with open("cube1_0-1.csv", "r") as g:
  raw_data = g.read()

# split the raw data to get a list of strings
list_o_strings = raw_data.split() 

# Use a nested list comprehension to step over the list of strings, splitting each string on its commas
# and converting each element (which will still be a string) to a floating point number
list_o_lists = [[float(num) for num in string.split(',')] for string in list_o_strings]

# make a numpy array
cube1 = np.array(list_o_lists)


# Function to make MatplotLib colormap
def gen_cmap(name, array, start, end):
    b3 = array[:,2] # value of blue at sample n
    b2 = array[:,2] # value of blue at sample n
    b1 = np.linspace(start, end, len(b2)) # position of sample n - ranges from 0 to 1
    
    # Setting up columns for tuples
    g3 = array[:,1]
    g2 = array[:,1]
    g1 = np.linspace(start, end, len(g2))
    
    r3 = array[:,0]
    r2 = array[:,0]
    r1 = np.linspace(start, end, len(r2))
    
    # Creating tuples
    R = sorted(zip(r1,r2,r3))
    G = sorted(zip(g1,g2,g3))
    B = sorted(zip(b1,b2,b3))
    
    # Transposing
    RGB = zip(R,G,B)
    rgb = zip(*RGB)
    
    # Creating dictionary
    k = ['red', 'green', 'blue']
    Cube1 = dict(zip(k,rgb))
    
    return clr.LinearSegmentedColormap(name, Cube1)
    
    
my_cmap = gen_cmap('my_colormap',cube1, 1, 0)
my_cmap_r = gen_cmap('my_colormap_r', np.flipud(cube1), 1, 0)


mycarta2=io.imread('https://mycarta.files.wordpress.com/2012/10/regio_distance_mineral_occurrences.png')

fig = plt.figure(figsize=(10, 8))

ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])

plt.imshow(mycarta2)
plt.show()


mycarta1=io.imread('https://mycarta.files.wordpress.com/2012/02/final-061.png')

fig = plt.figure(figsize=(10, 8))

ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])

plt.imshow(mycarta1)
plt.show()

data = np.loadtxt('continued.txt') # import Bouguer gravity residual map
print data.shape

vmin=data.min()
vmax=data.max()

fig = plt.figure(figsize=(10,7))

ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])

plt.imshow(data, cmap='gray', vmin=vmin, vmax=vmax)
plt.colorbar()
plt.show()

fig.savefig('data.png', dpi=200, bbox_inches='tight', pad_inches=0)


xc = np.loadtxt('xc.txt')
yc = np.loadtxt('yc.txt')

gauss = sp.ndimage.filters.gaussian_filter(data,1.) # apply mild smoothing

d = gauss.ravel() # rearrange data, X, and Y from 2D arrays to 1D arrays

y = xc.ravel() # x, y are switched because x in Fatiando is North-South. 
x = yc.ravel() # If we wanted to plot the data using coordinates we'd have to pass it, as , for example: contourf(y, x, ...)

xderiv = transform.derivx(x, y, d, data.shape) # calculate derivatives using Fatiando a Terra 
yderiv = transform.derivy(x, y, d, data.shape) 
zderiv = transform.derivz(x, y, d, data.shape)

xderiv2D = np.reshape(xderiv, (81, 81)) # reshape output arrays back to 2D
yderiv2D = np.reshape(yderiv, (81, 81))
zderiv2D = np.reshape(zderiv, (81, 81))

dx = xderiv2D # rename for convenience
dy = yderiv2D
dz = zderiv2D




fig = plt.figure(figsize=(10,7))

ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])

# gdz = sp.signal.medfilt2d(dz,3) # uncomment to apply very mild smoothing to remove a bit of station footprint

plt.imshow(dz, cmap='gray')
# plt.imshow(gdz, cmap='gray')

plt.colorbar()
plt.show()



fig = plt.figure(figsize=(20, 10)) # all three derivatives in one plot

ax3 = fig.add_subplot(1, 3, 1)
imshow(dz, cmap='gray')
ax3.set_xticks([])
ax3.set_yticks([])

ax4 = fig.add_subplot(1, 3, 2)
imshow(dx, cmap='gray')
ax4.set_xticks([])
ax4.set_yticks([])

ax5 = fig.add_subplot(1, 3, 3)
imshow(dy, cmap='gray')
ax5.set_xticks([])
ax5.set_yticks([])

tdx=np.sqrt(dx*dx + dy*dy) # calcualate total horizontal derivative

tdx01=(tdx-tdx.min())/(tdx.max()-tdx.min()) # rescale values to the range 0-1
tdxmin= tdx01.mean() - 2 * tdx01.std()      # calculate mean and standard deviation
tdxmax= tdx01.mean() + 2 * tdx01.std()      # to help wit hthe display of the data



fig = plt.figure(figsize=(12, 8))

ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])

plt.imshow(tdx01, cmap='bone_r', vmin = tdxmin, vmax = tdxmax) # display the total horizontal derivative
plt.colorbar()
plt.show()



tilt=np.arctan(dz/tdx) # calculate the tilt angle
gtilt = sp.ndimage.filters.gaussian_filter(tilt,1.) # apply minor gaussian smoothing

dxt,dyt = np.gradient(gtilt, 1, 1) # calculate THDR, total horizontal derivative of (smoothed) tilt angle
thdr=np.sqrt(dxt*dxt + dyt*dyt)

thdrmin= thdr.mean() - 2 * thdr.std()
thdrmax= thdr.mean() + 2 * thdr.std()

fig = plt.figure(figsize=(12, 8)) # display THDR

ax1 = fig.add_subplot(1, 1, 1)
ax1.set_xticks([])
ax1.set_yticks([])

plt.imshow(thdr,cmap='bone_r', vmin=thdrmin, vmax=thdrmax)
plt.colorbar()
plt.show()


theta=tdx/(np.sqrt(dx*dx + dy*dy +dz*dz)) # calculate theta map (e.g. cos(theta))
                                          # by dividing the total horizontal derivative by the analytical signal

thetamin= theta.mean() - 2 * theta.std()
thetamax= theta.mean() + 2 * theta.std()

# plot it

fig = plt.figure(figsize=(10,7))

ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])

plt.imshow(theta,cmap='bone_r', vmin=thetamin, vmax=thetamax)
plt.colorbar()
plt.show()

fig.savefig('theta.png', dpi=200)


nstd = np.loadtxt('nstd.txt') # import NSTD calculated using Matlab code from Cooper and Cowan.

nstdmin= nstd.mean() - 2 * nstd.std()
nstdmax= nstd.mean() + 2 * nstd.std()

fig = plt.figure(figsize=(10, 7)) # plot it 

ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])

plt.imshow(nstd, cmap='bone_r', vmin=nstdmin, vmax=nstdmax)
plt.colorbar()
plt.show()

z = (plt.contour(dz, [0.0])) # get the zero contours for dz
plt.close() # prevent Matplotlib from displaying the contours automatically (with default settings)

fig = plt.subplots(figsize=[5,5])
cmp = cm.get_cmap('bone', 2)
plt.contour(dz, [0.0], cmap=cmp, interpolation='none', origin="upper")  # plot with custom settings
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())


tdxn=np.real(np.arctan(tdx/np.absolute(dz))) # calculate normalized total horizontal derivative
tdxn = sp.signal.medfilt2d(tdxn,3) # apply very mild smoothing

tdxn01=(tdxn-tdxn.min())/(tdxn.max()-tdxn.min()) # rescale values to 0-1 range
tdxnmin= tdxn01.mean() - 2 * tdxn01.std()        # these help with display
tdxnmax= tdxn01.mean() + 2 * tdxn01.std()



fig = plt.figure(figsize=(12, 8))

ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])

plt.imshow(tdxn01,cmap='bone_r', vmin = tdxnmin, vmax=tdxnmax)
plt.colorbar()
plt.show()

fig = plt.figure(figsize=(12, 8))

ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])

plt.imshow(data,cmap=my_cmap_r)
plt.colorbar()
plt.show()

data_n=(data-data.min())/(data.max()-data.min())
img_array = plt.get_cmap(my_cmap_r)(data_n)

hsv = clr.rgb_to_hsv(img_array[:, :, :3])
hsv[:, :, 2] = 1-tdxn01
rgb = clr.hsv_to_rgb(hsv)


fig = plt.figure(figsize=(12, 8))

ax1 = fig.add_subplot(1, 1, 1)
ax1.set_xticks([])
ax1.set_yticks([])

plt.imshow(rgb)
plt.show()

line = z.collections[0].get_paths()[2].vertices  # selected 2nd by trial and error to get the longest contour segment

fig = plt.figure(figsize=(7, 5))

ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])

plt.imshow(data,cmap=my_cmap_r, origin="upper")
plt.plot(line[:,0], line[:,1], 'k-')

plt.colorbar()
plt.show()

fig = plt.figure(figsize=(12, 8))

ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])
plt.imshow(data,cmap=my_cmap_r, origin="upper")

for collection in z.collections:
    paths = collection.get_paths()
    for path in paths:
        line=(path.vertices)
        plt.plot(line[:,0], line[:,1], 'k-')
plt.colorbar()
plt.show()

thresh_th = threshold_otsu(theta)
binary_th = theta > thresh_th

# normalized total horizontal derivative
thresh_td = threshold_otsu(tdxn01)
binary_td = tdxn01 > thresh_td

# nstd
thresh_ns = threshold_otsu(nstd)
binary_ns = nstd > thresh_ns


# plot them 
fig = plt.figure(figsize=(20, 7))

ax1 = fig.add_subplot(1, 3, 1)
ax1.set_xticks([])
ax1.set_yticks([])
imshow(binary_th, cmap=plt.cm.gray)

ax2 = fig.add_subplot(1, 3, 2)
ax2.set_xticks([])
ax2.set_yticks([])
imshow(binary_td, cmap=plt.cm.gray)

ax3 = fig.add_subplot(1, 3, 3)
ax3.set_xticks([])
ax3.set_yticks([])
imshow(binary_ns, cmap=plt.cm.gray)


binary_th = color.rgb2gray(theta) > 0.9

# normalized total horizontal derivative
binary_td = color.rgb2gray(tdxn01) > 0.75

# nstd
binary_ns = color.rgb2gray(nstd) > 0.45


# plot them 
fig = plt.figure(figsize=(20, 7))

ax1 = fig.add_subplot(1, 3, 1)
ax1.set_xticks([])
ax1.set_yticks([])
imshow(binary_th, cmap=plt.cm.gray)

ax2 = fig.add_subplot(1, 3, 2)
ax2.set_xticks([])
ax2.set_yticks([])
imshow(binary_td, cmap=plt.cm.gray)

ax3 = fig.add_subplot(1, 3, 3)
ax3.set_xticks([])
ax3.set_yticks([])
imshow(binary_ns, cmap=plt.cm.gray)


label_objects_th, nb_labels_th = sp.ndimage.label(binary_th) # label all white objects ( made up of ones)
sizes_th = np.bincount(label_objects_th.ravel())             # calculate every labeled object's size
mask_sizes_th = sizes_th > 75                                # mask with zeros all objects smaller than
mask_sizes_th[0] = 0                                         # the given size
binary_cleaned_th = mask_sizes_th[label_objects_th]

# normalized total horizontal derivative
label_objects_td, nb_labels_td = sp.ndimage.label(binary_td)
sizes_td = np.bincount(label_objects_td.ravel())
mask_sizes_td = sizes_td > 75
mask_sizes_td[0] = 0
binary_cleaned_td = mask_sizes_td[label_objects_td]

# nstd
label_objects_ns, nb_labels_ns = sp.ndimage.label(binary_ns)
sizes_ns = np.bincount(label_objects_ns.ravel())
mask_sizes_ns = sizes_ns > 75
mask_sizes_ns[0] = 0
binary_cleaned_ns = mask_sizes_ns[label_objects_ns]


# plot 
fig = plt.figure(figsize=(20, 7))

ax1 = fig.add_subplot(1, 3, 1)
ax1.set_xticks([])
ax1.set_yticks([])
imshow(binary_cleaned_th, cmap=plt.cm.gray)

ax2 = fig.add_subplot(1, 3, 2)
ax2.set_xticks([])
ax2.set_yticks([])
imshow(binary_cleaned_td, cmap=plt.cm.gray)

ax3 = fig.add_subplot(1, 3, 3)
ax3.set_xticks([])
ax3.set_yticks([])
imshow(binary_cleaned_ns, cmap=plt.cm.gray)


selem = disk(1)

closed_th = closing(binary_cleaned_th, selem)
closed_td = closing(binary_cleaned_td, selem)
closed_ns = closing(binary_cleaned_ns, selem)

eroded_th = erosion(closed_th, selem)
eroded_td = erosion(closed_td, selem)
eroded_ns = erosion(closed_ns, selem)

# plot them 
fig = plt.figure(figsize=(20, 7))

ax1 = fig.add_subplot(1, 3, 1)
ax1.set_xticks([])
ax1.set_yticks([])
imshow(closed_th, cmap=plt.cm.gray)

ax2 = fig.add_subplot(1, 3, 2)
ax2.set_xticks([])
ax2.set_yticks([])
imshow(closed_td, cmap=plt.cm.gray)

ax3 = fig.add_subplot(1, 3, 3)
ax3.set_xticks([])
ax3.set_yticks([])
imshow(closed_ns, cmap=plt.cm.gray)






skeleton_th=skeletonize(binary_th)
skeleton_cleaned_th=skeletonize(binary_cleaned_th)
skeleton_cleaned_th1=skeletonize(closed_th)

skeleton_td=skeletonize(binary_td)
skeleton_cleaned_td=skeletonize(binary_cleaned_td)
skeleton_cleaned_td1=skeletonize(closed_td)

skeleton_ns=skeletonize(binary_ns)
skeleton_cleaned_ns=skeletonize(binary_cleaned_ns)
skeleton_cleaned_ns1=skeletonize(closed_ns)
skeleton_cleaned_ns1=np.pad(skeleton_cleaned_ns1, 4, 'minimum')  # padding the Matlab NSTD array since it is smaller


fig = plt.figure(figsize=(20, 7))

ax4 = fig.add_subplot(1, 3, 1)
imshow(skeleton_td, cmap='bone_r', interpolation='none')

ax5 = fig.add_subplot(1, 3, 2)
imshow(skeleton_cleaned_td, cmap='bone_r', interpolation='none')

ax6 = fig.add_subplot(1, 3, 3)
imshow(skeleton_cleaned_td1, cmap='bone_r', interpolation='none')

ax4.set_xticks([])
ax4.set_yticks([])

ax5.set_xticks([])
ax5.set_yticks([])

ax6.set_xticks([])
ax6.set_yticks([])



fig = plt.figure(figsize=(20, 7))

ax4 = fig.add_subplot(1, 3, 1)
imshow(skeleton_cleaned_td1, cmap='bone_r', interpolation='none')

ax5 = fig.add_subplot(1, 3, 2)
imshow(skeleton_cleaned_th1, cmap='bone_r', interpolation='none')

ax6 = fig.add_subplot(1, 3, 3)
imshow(skeleton_cleaned_ns1, cmap='bone_r', interpolation='none')

ax4.set_xticks([])
ax4.set_yticks([])

ax5.set_xticks([])
ax5.set_yticks([])

ax6.set_xticks([])
ax6.set_yticks([])



w = 4.05 # define figure size as a 20th of the size of the array
h = 4.05

cmp = cm.get_cmap('bone', 2) # grab a two-color colormap

fig = plt.figure(frameon=False)
fig.set_size_inches(w,h)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.set_xticks([])
ax.set_yticks([])

plt.imshow(thdr,cmap='bone_r', vmin=thdrmin, vmax=thdrmax)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
fig.savefig('THDR.png', dpi=20) # export at 20 dpi to get an image of exactly the same size of the array

inTHDR=io.imread('THDR.png')       # import back THDR image
binary_thdr = ~(color.rgb2gray(inTHDR) > 0.34)*1.0  # binarize it; the result of rgb2gray is logical (True, False)
                                                    # hence the multiplication by 1.0
fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])
plt.imshow(binary_thdr, cmap='bone')
plt.show()


label_objects_thdr1, nb_labels_thdr1 = sp.ndimage.label(binary_thdr) # clean up small objects thdr
sizes_thdr1 = np.bincount(label_objects_thdr1.ravel())
mask_sizes_thdr1 = sizes_thdr1 > 15
mask_sizes_thdr1[0] = 0
binary_cleaned_thdr1 = mask_sizes_thdr1[label_objects_thdr1]

fig = plt.figure(figsize=(20, 7))
ax1 = fig.add_subplot(1, 3, 1)
ax1.set_xticks([])
ax1.set_yticks([])
imshow(binary_cleaned_thdr1, cmap='bone')


closed_thdr = closing(binary_cleaned_thdr1, selem)

# plot it 
fig = plt.figure(figsize=(5, 5))

ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])

plt.imshow(closed_thdr, cmap='bone')
plt.show()

skeleton_cleaned_THDR=skeletonize(closed_thdr)

fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(1, 1, 1)
imshow(skeleton_cleaned_THDR, cmap='bone_r',interpolation='none')
ax.set_xticks([])
ax.set_yticks([])


selem2 = disk(1)
dilated_skeleton_td1 = dilation (skeleton_cleaned_td1, selem2)
dilated_skeleton_thdr1 = dilation (skeleton_cleaned_THDR, selem2)
dilated_skeleton_ns1 = dilation (skeleton_cleaned_ns1, selem2)


dilated_skeleton_td1 = dilated_skeleton_td1*1.0
dilated_skeleton_thdr1 = dilated_skeleton_thdr1*1.0
dilated_skeleton_ns1 = dilated_skeleton_ns1*1.0

print np.amax(dilated_skeleton_td1), np.amin(dilated_skeleton_td1)
print np.amax(dilated_skeleton_thdr1), np.amin(dilated_skeleton_thdr1)
print np.amax(dilated_skeleton_ns1), np.amin(dilated_skeleton_ns1)


fig = plt.figure(figsize=(20, 8))

ax1 = fig.add_subplot(1, 3, 1)
imshow(dilated_skeleton_td1, cmap='bone_r',interpolation='spline16')
ax1.set_xticks([])
ax1.set_yticks([])

ax2 = fig.add_subplot(1, 3, 2)
imshow(dilated_skeleton_thdr1, cmap='bone_r',interpolation='spline16')
ax2.set_xticks([])
ax2.set_yticks([])

ax3 = fig.add_subplot(1, 3, 3)
imshow(dilated_skeleton_ns1, cmap='bone_r',interpolation='spline16')
ax3.set_xticks([])
ax3.set_yticks([])


collocation =  dilated_skeleton_thdr1 + dilated_skeleton_td1 + dilated_skeleton_ns1


ch = plt.get_cmap('afmhot')

# extract the colormap RGBA values at the 16 points
rgba = ch(np.arange(256))

# slice rgba to discard Alpha
rgb01= rgba[0:180,:3]

# setting up color arrays
r1 = rgb01[:, 0] # value of Red for the nth sample
g1 = rgb01[:, 1] # value of Green for the nth sample
b1 = rgb01[:, 2] # value of Blue for the nth sample

r2 = r1 # value of Red at the nth sample
r0 = np.linspace(0, 1, len(r1)) # position of the nth Red sample within the range 0 to 1

g2 = g1 # value of Green at the nth sample
g0 = np.linspace(0, 1, len(g1)) # position of the nth Green sample within the range 0 to 1

b2 = b1 # value of Blue at the nth sample
b0 = np.linspace(0, 1, len(b1)) # position of the nth Blue sample within the range 0 to 1

# creating lists
R = zip(r0, r1, r2)
G = zip(g0, g1, g2)
B = zip(b0, b1, b2)

# creating list of above lists and transposing
RGB = zip(R, G, B)
rgb = zip(*RGB)
#print rgb

# creating dictionary
k=['red', 'green', 'blue'] # makes list of keys
my_hot=dict(zip(k,rgb)) # makes a dictionary from list of keys and list of values
m_hot = clr.LinearSegmentedColormap('my_hot', my_hot)


fig = plt.figure(figsize=(10, 7))

ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])

plt.imshow(collocation,cmap=m_hot, interpolation='spline16')
cb=plt.colorbar()
cb.set_ticks([]) # taking off the labels from the colormap to add custom ones in image editor

fig.savefig('collocation_full.png', dpi=200, bbox_inches='tight', pad_inches=0)


finalmap = np.ma.masked_where(collocation <2., collocation)
finalmap = finalmap.filled(0)

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])

plt.imshow(finalmap,cmap=m_hot, vmax=3, vmin=0, interpolation='spline16')
cb = plt.colorbar()
cb.set_ticks([])

fig.savefig('colocation.png', dpi=200, bbox_inches='tight', pad_inches=0)


fig = plt.figure(figsize=(10, 7))

ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])

plt.imshow(tdxn,  cmap='Greys')
finalmap[finalmap==0] = np.nan
plt.imshow(finalmap, cmap='afmhot', vmin=0, vmax=4)
plt.show()


fig.savefig('tdxn_and_colocation.png', dpi=200, bbox_inches='tight', pad_inches=0)


Fig3=io.imread('Figure3_big.png')

fig = plt.figure(figsize=(20, 14))
ax = fig.add_subplot(1, 1, 1)
ax.set_xticks([])
ax.set_yticks([])
plt.imshow(Fig3)
plt.show()







































































    
    
    
