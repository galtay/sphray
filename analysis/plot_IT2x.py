import sphray_io
import numpy as np

import matplotlib as mpl

mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['font.size'] = 20.0
mpl.rcParams['xtick.labelsize'] = 20.0
mpl.rcParams['ytick.labelsize'] = 20.0

import matplotlib.pyplot as plt

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

#====================================================================
# define input parameters
#====================================================================

snapdir  = "../../sphray_output/IT2_N64/r6"
snapbase = "snap"
snapnum = 1
snapnumstr = '{0:03}'.format( snapnum )

sfile  = snapdir + "/" + snapbase + "_" + snapnumstr
pngfile = "T2x" + snapnumstr + ".png"
cmpfile= "CmpData/CmpT2_" + snapnumstr + "x.txt"

print
print "comparison project file:", cmpfile
print "sphray output file:", sfile
print "png file if doing file output:", pngfile
print


#====================================================================
# define some physical constants
#====================================================================

PROTONMASS = 1.6726e-24 # proton mass [g]
GLEN  = 3.085678e21     # [cm h^-1]
GMASS = 1.989e43        # [g h^-1]
GVEL  = 1.0e5           # [cm s^-1]

GTIME = GLEN / GVEL                 # [s h^-1]
GRHO  = GMASS / GLEN**3             # [g cm^-3 h^2]
GPRS  = GMASS / GLEN / GTIME**2     # [g cm^-1 s^-2 h^2]
GENRG = GMASS * GLEN**2 / GTIME**2  # [erg h^-1]
GLUM  = GENRG / GTIME               # [erg/s]


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


#====================================================================
# read comparison project data
#====================================================================

codes = ['xx', 'c2ray', 'otvet', 'crash', 'rsph', 'art', 'ftte', 
         'zeus', 'ift']

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

print 'finished reading and sorting' 

fig = plt.figure( figsize=(10,10) ) 
ax = fig.add_subplot(111)

ax.scatter( sdata['rad']/(shead['boxlen'][0]/2), 
             np.log10( sdata['xHI'] ), s=30, 
            facecolors='grey', edgecolors='grey', alpha=0.1)

ax.scatter( sdata['rad']/(shead['boxlen'][0]/2), 
             np.log10( 1.0-sdata['xHI'] ), s=30, 
            facecolors='grey', edgecolors='grey', alpha=0.1)

for i,c in enumerate(codes[1:]):

    this_line = next(linecycler)
    color = this_line['color']
    ls = this_line['ls']


    ax.plot( cdata['xx'], np.log10( cdata[c] ), label=c, 
             lw=2.0, color=color, ls=ls )


    ax.plot( cdata['xx'], np.log10( 1.0-cdata[c] ),  
             lw=2.0, color=color, ls=ls )


ax.legend(loc='lower center', ncol=2)

ax.set_xlabel( r'$r/L_{\rm box}$', fontsize=20 )
ax.set_xlim( -0.05, 1.05 )

ax.set_ylabel( r'$x_{\rm HI}$', fontsize=20 )
ax.set_ylim( -5.3, 0.1 )

fig.savefig( pngfile )


