import sphray_io
import numpy as np
import matplotlib.pyplot as plt


#====================================================================
# define input parameters
#====================================================================

snapdir  = "../../sphray_output/IT1_N64/r6"
snapbase = "snap"
snapnum = 1
snapnumstr = '{0:03}'.format( snapnum )

sfile  = snapdir + "/" + snapbase + "_" + snapnumstr
pngfile = "T1x" + snapnumstr + ".png"
cmpfile= "CmpData/CmpT1_" + snapnumstr + "x.txt"

print
print "comparison project file:", cmpfile
print "sphray output file:", sfile
print "png file if doing file output:", pngfile
print




#====================================================================
# read sphray snapshot
#====================================================================

# create a file object and read data
#---------------------------------------------------

sf = sphray_io.SphrayFile()

shead = sf.read_header(sfile)
sdata = sf.read_data_1(sfile) 

# sort particles by radius
#---------------------------------------------------

sdata['pos'] = sdata['pos'] - shead['boxlen'][0]/2
sdata['rad'] = np.sqrt( sdata['pos'][:,0]**2 +  
                        sdata['pos'][:,1]**2 +  
                        sdata['pos'][:,2]**2 )

print 'min/max pos: ', sdata['pos'].min(), sdata['pos'].max()
print 'min/max rad: ', sdata['rad'].min(), sdata['rad'].max()

sdata = sf.convert_data_to_structured_array(sdata)

sdata.sort(order='rad')


# bin particles radially with an equal number of 
# particles per bin
#---------------------------------------------------

nbins = 1000
ngas = sdata['rad'].size

nperbin = ngas / nbins

class SphRadialAverage:
    pass

rav = SphRadialAverage()

rav.xx = np.zeros( nbins )
rav.xHI_mean = np.zeros( nbins )
rav.xHI_median = np.zeros( nbins )


for i in xrange(nbins):
    ii = i * nperbin
    ff = ii + nperbin
    if i==nbins-1:
        ff = ngas

    rads = sdata['rad'][ii:ff]
    xHI = sdata['xHI'][ii:ff]

    rav.xx[i] = np.mean(rads)
    rav.xHI_mean[i] = np.mean( np.log10(xHI) )
    rav.xHI_median[i] = np.median( np.log10(xHI) )


#====================================================================
# read comparison project data
#====================================================================

codes = ['xx', 'c2ray', 'otvet', 'crash', 'rsph', 'art', 'ftte', 
         'simplex', 'zeus', 'flash', 'ift']

Ncmpbins = 100
cdata_pre = np.loadtxt( cmpfile )


cdata={}
for i,c in enumerate(codes): 
    ii = i * Ncmpbins
    ff = ii + 100
    dat = cdata_pre[ii:ff]
    cdata[c] = dat



#====================================================================
# make plot
#====================================================================



# define linestyles 
#--------------------------------------

from itertools import cycle

lines = [
    {'color':'blue', 'ls':'-'}, 
    {'color':'green', 'ls':'-'}, 
    {'color':'red', 'ls':'-'}, 
    {'color':'gold', 'ls':'-'}, 
    {'color':'cyan', 'ls':'-'}, 
    {'color':'purple', 'ls':'-'}, 
    {'color':'black', 'ls':'-'}, 
    {'color':'red', 'ls':'--'}, 
    {'color':'green', 'ls':'--'}, 
    
]
linecycler = cycle(lines)


# set up figure
#--------------------------------------

fig = plt.figure( figsize=(10,10) ) 
ax = fig.add_subplot(111)

# plot SPHRAY results
#--------------------------------------

ax.scatter( sdata['rad']/(shead['boxlen'][0]/2), 
             np.log10( sdata['xHI'] ), s=30, 
            facecolors='grey', edgecolors='grey', alpha=0.1)

ax.scatter( sdata['rad']/(shead['boxlen'][0]/2), 
             np.log10( 1.0-sdata['xHI'] ), s=30, 
            facecolors='grey', edgecolors='grey', alpha=0.1)

ax.plot( rav.xx / (shead['boxlen'][0]/2), 
         rav.xHI_median, color='lime', lw=3.0,
         alpha=0.7, zorder=4)

ax.plot( rav.xx / (shead['boxlen'][0]/2), 
         np.log10( 1.0 - 10**rav.xHI_median), 
         color='lime', lw=3.0,
         alpha=0.7, zorder=4)



# plot comparison project results
#--------------------------------------


for i,c in enumerate(codes[1:]):

    this_line = next(linecycler)
    color = this_line['color']
    ls = this_line['ls']


    ax.plot( cdata['xx'], np.log10( cdata[c] ), label=c, 
             lw=2.0, color=color, ls=ls )


    ax.plot( cdata['xx'], np.log10( 1.0-cdata[c] ),  
             lw=2.0, color=color, ls=ls )



# finalize plot
#--------------------------------------


ax.legend(loc='lower center', ncol=2)

ax.set_xlabel( r'$r/L_{\rm box}$', fontsize=20 )
ax.set_xlim( -0.05, 1.05 )

ax.set_ylabel( r'$x_{\rm HI}$', fontsize=20 )
ax.set_ylim( -5.3, 0.1 )

fig.savefig( pngfile )



