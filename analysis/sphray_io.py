import numpy as np


class SphrayFile:

    """ A class to handle IO with SPHRAY output files. """

    header_dtypes = np.dtype([
            ('head','i4',1),
            ('npar_file','i4',6),
            ('mass_arr','f8', 6),
            ('time', 'f8',1),
            ('z', 'f8', 1),
            ('flag_sfr', 'i4', 1),
            ('flag_feedback', 'i4', 1),
            ('npar_snap', 'i4', 6),
            ('flag_cool', 'i4', 1),
            ('nfiles', 'i4', 1),
            ('boxlen', 'f8', 1),
            ('omegaM', 'f8', 1),
            ('omegaL', 'f8', 1),
            ('h', 'f8', 1),
            ('flag_age', 'i4', 1),
            ('flag_metals', 'i4', 1),
            ('npar_hw', 'i4', 6),
            ('flag_entr_ics', 'i4', 1),
            ('omegaB', 'f8', 1),
            ('nrays', 'i8', 1),
            ('flag_Hmf', 'i4', 1),
            ('flag_Hemf', 'i4', 1),
            ('flag_helium', 'i4', 1),
            ('flag_gammaHI', 'i4', 1),
            ('flag_cloudy', 'i4', 1),
            ('flag_eos', 'i4', 1),
            ('unused', 'i4', 5),
            ('tail', 'i4', 1),
            ])



    data_keys = ['pos', 'vel', 'ids', 'mass', 'u', 'rho', 'ye',
                 'xHI', 'hsml', 'T', 'Hmf', 'Hemf', 'xHeI', 'xHeII', 
                 'GHI', 'xHI_cloudy', 'eos', 'lasthit']




    data_type_particle = {
        'pos': ( 'f4', (3) ), 
        'vel': ( 'f4', (3) ), 
        'ids': ( 'i4', (1) ), 
        'mass': ( 'f4', (1) ), 
        'u': ( 'f4', (1) ), 
        'rho': ( 'f4', (1) ), 
        'ye': ( 'f4', (1) ),
        'xHI': ( 'f4', (1) ), 
        'hsml': ( 'f4', (1) ), 
        'T': ( 'f4', (1) ), 
        'Hmf': ( 'f4', (1) ), 
        'Hemf': ( 'f4', (1) ), 
        'xHeI': ( 'f4', (1) ), 
        'xHeII': ( 'f4', (1) ), 
        'GHI': ( 'f4', (1) ), 
        'xHI_cloudy': ( 'f4', (1) ), 
        'eos': ( 'f4', (1) ), 
        'lasthit': ( 'i8', (1) )
        }
    



    def read_data_block( self, f, key ):
        head = np.fromfile( f, dtype='i4', count=1 )
        tmp = np.fromfile( f, dtype=self.data_types[key], count=1 ) 
        tail = np.fromfile( f, dtype='i4', count=1 )            
        assert( head[0]==tail[0] )
        shape = tmp.shape[1:]
        tmp = tmp.reshape( shape )
        return tmp



    def set_data_types( self, ngas ):

        self.data_types = {'pos': np.dtype( ('f4', (ngas,3)) ), 
                           'vel': np.dtype( ('f4', (ngas,3)) ), 
                           'ids': np.dtype( ('i4', (ngas)) ), 
                           'mass': np.dtype( ('f4', (ngas)) ), 
                           'u': np.dtype( ('f4', (ngas)) ), 
                           'rho': np.dtype( ('f4', (ngas)) ), 
                           'ye': np.dtype( ('f4', (ngas)) ),
                           'xHI': np.dtype( ('f4', (ngas)) ), 
                           'hsml': np.dtype( ('f4', (ngas)) ), 
                           'T': np.dtype( ('f4', (ngas)) ), 
                           'Hmf': np.dtype( ('f4', (ngas)) ), 
                           'Hemf': np.dtype( ('f4', (ngas)) ), 
                           'xHeI': np.dtype( ('f4', (ngas)) ), 
                           'xHeII': np.dtype( ('f4', (ngas)) ), 
                           'GHI': np.dtype( ('f4', (ngas)) ), 
                           'xHI_cloudy': np.dtype( ('f4', (ngas)) ), 
                           'eos': np.dtype( ('f4', (ngas)) ), 
                           'lasthit': np.dtype( ('i8', (ngas)) )}
        

        



    def set_data_exists( self, header ):

        self.data_exists = {'pos': True, 'vel': True, 'ids': True, 
                            'mass': True, 'u': True, 'rho': True, 
                            'ye': True, 'xHI': True, 'hsml': True, 
                            'T': True, 'Hmf': False, 'Hemf': False, 
                            'xHeI': False, 'xHeII': False, 'GHI': False, 
                            'xHI_cloudy': False, 'eos': False, 'lasthit': True} 
        
        if header['flag_Hmf'][0]:
            self.data_exists['Hmf'] = True

        if header['flag_Hemf'][0]:
            self.data_exists['Hemf'] = True

        if header['flag_helium'][0]:
            self.data_exists['xHeI'] = True
            self.data_exists['xHeII'] = True

        if header['flag_gammaHI'][0]:
            self.data_exists['GHI'] = True

        if header['flag_cloudy'][0]:
            self.data_exists['xHI_cloudy'] = True

        if header['flag_eos'][0]:
            self.data_exists['eos'] = True






    def read_header( self, fname ):
        header=None
        with open(fname, 'rb') as f:
            header = np.fromfile( f, dtype=self.header_dtypes, count=1 )
            assert( header['head']==header['tail'] )
        return header



    def read_data_1( self, fname ):

        data=None

        with open(fname, 'rb') as f:
            header = np.fromfile( f, dtype=self.header_dtypes, count=1 )
            assert( header['head']==header['tail'] )
            ngas = header['npar_file'][0][0]
            print 'ngas = ', ngas

            self.set_data_exists( header )
            self.set_data_types( ngas )

            data = {}
            for key in self.data_keys:
                
                if self.data_exists[key]:
                    tmp = self.read_data_block( f, key ) 
                    data[key] = tmp

        return data



    def convert_data_to_structured_array( self, data ):
        """ Converts the dictionary representation of data (i.e. key = 'pos'
        leads to all values of position) into a structured array version of
        data (i.e. an array of particle records that can be accessed by 
        particle number or by data tag). """
        
        dtype = []
        for key in data.keys():            
            if key in self.data_type_particle.keys():
                dtype.append( (key, 
                               self.data_type_particle[key][0],
                               self.data_type_particle[key][1]) )
            else:
                dtype.append( (key,'f4',(1)) )

        npar = data[data.keys()[0]].shape[0]
        sa = np.zeros( npar, dtype=dtype )

        for key in data.keys():
            sa[key] = data[key]

        return sa
