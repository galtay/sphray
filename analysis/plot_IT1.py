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

print 'finished reading and sorting' 

fig = plt.figure() 
ax = fig.add_subplot(111)

ax.scatter( sdata['rad']/(shead['boxlen'][0]/2), 
             np.log10( sdata['xHI'] ), s=10 )

for i,c in enumerate(codes[1:]):
    ax.plot( cdata['xx'], np.log10( cdata[c] ) )


ax.set_xlim( -0.05, 1.05 )
ax.set_ylim( -5.3, 0.1 )




