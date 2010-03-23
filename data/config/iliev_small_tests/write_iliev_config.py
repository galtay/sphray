config_list= ['Verbosity', 'DoTestScenario', 'TestScenario', 'JustInit',
              'Comoving', 'IsoTemp', 'FixSnapTemp', 'EOStemp', 'InitxHI',
              'RayDepletion', 'IntSeed', 'StaticFieldSimTime',
              'StaticSimTimeUnit', 'InputType', 'SnapPath', 'SourcePath',
              'SpectraFile', 'b2cdFile', 'AtomicRatesFile', 'ParFileBase',
              'SourceFileBase', 'StartSnapNum', 'EndSnapNum',
              'ParFilesPerSnap', 'SourceFilesPerSnap', 'RayScheme',
              'ForcedRayNumber', 'RayStats', 'BndryCond', 'RecRayTol',
              'RayPhotonTol', 'WallSampling', 'OnTheSpotH', 'OnTheSpotHe',
              'HydrogenCaseA', 'HeliumCaseA', 'IonTempSolver', 'Tfloor',
              'Tceiling', 'xfloor', 'xceiling', 'NeBackground',
              'NraysUpdateNoHits', 'RecRaysPerSrcRay', 'H_mf', 'He_mf',
              'OutputDir', 'OutputFileBase', 'OutputType', 'OutputTiming',
              'NumStdOuts', 'DoInitialOutput', 'IonFracOutRays',
              'ForcedOutFile', 'ForcedUnits', 'PartPerCell']




defaults= ['3', 'true', 'iliev_test1', 'False',
           'False', '1.0d4', 'False', '1.0d3', '-99.0d0',
           'True', '0123456789', '500.0d0',
           'myr', '1', '../../rtcp_snapshots', '../../rtcp_snapshots',
           '../data/spectra/thermal1e5.cdf',
           '../data/column_depth/latmon_b2cd_table.txt',
           '../data/atomic_rates/atomic_rates_Hui.txt',
           'gadget_glass_L13.2_N64', 'test1_sources',
           '001', '001', '1', '1', 'raynum', '1.0d6', 'False', '0', '3.0d-1',
           '1.0d-10', '3', 'True', 'True', 'False', 'False', '1',
           '1.0d0', '1.0d9', '0.0d0', '1.0d0', '1.0d-8',
           '0','10', '1.0d0', '0.0d0',
           '../../sphray_output/IT1/r6', 'snap', '1', 'forced',
           '0', 'False', '10000',
           '../../rtcp_snapshots/test1_output_times.txt', 'myr', '12']





config_files= ['iliev_test1_N64_R6.config',           'iliev_test1_N64_R7.config',           'iliev_test1_N64_R8.config',
               'iliev_test1_N128_R6.config',          'iliev_test1_N128_R7.config',          'iliev_test1_N128_R8.config',
               'iliev_test1_He_N64_R6.config',        'iliev_test1_He_N64_R7.config',        'iliev_test1_He_N64_R8.config',
               'iliev_test1_He_N128_R6.config',       'iliev_test1_He_N128_R7.config',       'iliev_test1_He_N128_R8.config',
               'iliev_test1_HM01QG_N64_R6.config',    'iliev_test1_HM01QG_N64_R7.config',    'iliev_test1_HM01QG_N64_R8.config',
               'iliev_test1_HM01QG_N128_R6.config',   'iliev_test1_HM01QG_N128_R7.config',   'iliev_test1_HM01QG_N128_R8.config',
               'iliev_test1_HM01QGnd_N64_R6.config',  'iliev_test1_HM01QGnd_N64_R7.config',  'iliev_test1_HM01QGnd_N64_R8.config',
               'iliev_test1_HM01QGnd_N128_R6.config', 'iliev_test1_HM01QGnd_N128_R7.config', 'iliev_test1_HM01QGnd_N128_R8.config',
               'iliev_test2_N64_R6.config',           'iliev_test2_N64_R7.config',           'iliev_test2_N64_R8.config',
               'iliev_test2_N128_R6.config',          'iliev_test2_N128_R7.config',          'iliev_test2_N128_R8.config']

dotests= ['True', 'True', 'True',
          'True', 'True', 'True',
          'True', 'True', 'True',
          'True', 'True', 'True',
          'False', 'False', 'False',
          'False', 'False', 'False',
          'False', 'False', 'False',
          'False', 'False', 'False',
          'True', 'True', 'True',
          'True', 'True', 'True']

test_names= ['iliev_test1', 'iliev_test1', 'iliev_test1',
             'iliev_test1', 'iliev_test1', 'iliev_test1',
             'iliev_test1He', 'iliev_test1He', 'iliev_test1He',
             'iliev_test1He', 'iliev_test1He', 'iliev_test1He',
             'none', 'none', 'none',
             'none', 'none', 'none',
             'none', 'none', 'none',
             'none', 'none', 'none',
             'iliev_test2', 'iliev_test2', 'iliev_test2',
             'iliev_test2', 'iliev_test2', 'iliev_test2']

iso_temps= ['1.0d4', '1.0d4', '1.0d4',
            '1.0d4', '1.0d4', '1.0d4',
            '1.0d4', '1.0d4', '1.0d4',
            '1.0d4', '1.0d4', '1.0d4',
            '1.0d4', '1.0d4', '1.0d4',
            '1.0d4', '1.0d4', '1.0d4',
            '1.0d4', '1.0d4', '1.0d4',
            '1.0d4', '1.0d4', '1.0d4',
            '0.0d0', '0.0d0', '0.0d0',
            '0.0d0', '0.0d0', '0.0d0']

deplete_rays= ['True', 'True', 'True',
               'True', 'True', 'True',
               'True', 'True', 'True',
               'True', 'True', 'True',
               'True', 'True', 'True',
               'True', 'True', 'True',
               'False', 'False', 'False',
               'False', 'False', 'False',
               'True', 'True', 'True',
               'True', 'True', 'True']

spectra_files= ['../data/spectra/thermal1e5.cdf',        '../data/spectra/thermal1e5.cdf',        '../data/spectra/thermal1e5.cdf',
                '../data/spectra/thermal1e5.cdf',        '../data/spectra/thermal1e5.cdf',        '../data/spectra/thermal1e5.cdf',
                '../data/spectra/thermal1e5.cdf',        '../data/spectra/thermal1e5.cdf',        '../data/spectra/thermal1e5.cdf',
                '../data/spectra/thermal1e5.cdf',        '../data/spectra/thermal1e5.cdf',        '../data/spectra/thermal1e5.cdf',
                '../data/spectra/hm01/hm01qg_z0.00.cdf', '../data/spectra/hm01/hm01qg_z0.00.cdf', '../data/spectra/hm01/hm01qg_z0.00.cdf',
                '../data/spectra/hm01/hm01qg_z0.00.cdf', '../data/spectra/hm01/hm01qg_z0.00.cdf', '../data/spectra/hm01/hm01qg_z0.00.cdf',
                '../data/spectra/hm01/hm01qg_z0.00.cdf', '../data/spectra/hm01/hm01qg_z0.00.cdf', '../data/spectra/hm01/hm01qg_z0.00.cdf',
                '../data/spectra/hm01/hm01qg_z0.00.cdf', '../data/spectra/hm01/hm01qg_z0.00.cdf', '../data/spectra/hm01/hm01qg_z0.00.cdf',
                '../data/spectra/thermal1e5.cdf',        '../data/spectra/thermal1e5.cdf',        '../data/spectra/thermal1e5.cdf',
                '../data/spectra/thermal1e5.cdf',        '../data/spectra/thermal1e5.cdf',        '../data/spectra/thermal1e5.cdf']

parbases= ['gadget_glass_L13.2_N64',  'gadget_glass_L13.2_N64',  'gadget_glass_L13.2_N64',
           'gadget_glass_L13.2_N128', 'gadget_glass_L13.2_N128', 'gadget_glass_L13.2_N128',
           'gadget_glass_L13.2_N64',  'gadget_glass_L13.2_N64',  'gadget_glass_L13.2_N64',
           'gadget_glass_L13.2_N128', 'gadget_glass_L13.2_N128', 'gadget_glass_L13.2_N128',
           'gadget_glass_L13.2_N64',  'gadget_glass_L13.2_N64',  'gadget_glass_L13.2_N64',
           'gadget_glass_L13.2_N128', 'gadget_glass_L13.2_N128', 'gadget_glass_L13.2_N128',
           'gadget_glass_L13.2_N64',  'gadget_glass_L13.2_N64',  'gadget_glass_L13.2_N64',
           'gadget_glass_L13.2_N128', 'gadget_glass_L13.2_N128', 'gadget_glass_L13.2_N128',
           'gadget_glass_L13.2_N64',  'gadget_glass_L13.2_N64',  'gadget_glass_L13.2_N64',
           'gadget_glass_L13.2_N128', 'gadget_glass_L13.2_N128', 'gadget_glass_L13.2_N128']

srcbases= ['test1_sources', 'test1_sources', 'test1_sources',
           'test1_sources', 'test1_sources', 'test1_sources',
           'test2_sources', 'test2_sources', 'test2_sources',
           'test2_sources', 'test2_sources', 'test2_sources',
           'test1_HM01QG_sources', 'test1_HM01QG_sources', 'test1_HM01QG_sources',
           'test1_HM01QG_sources', 'test1_HM01QG_sources', 'test1_HM01QG_sources',
           'test1_HM01QG_sources', 'test1_HM01QG_sources', 'test1_HM01QG_sources',
           'test1_HM01QG_sources', 'test1_HM01QG_sources', 'test1_HM01QG_sources',
           'test2_sources', 'test2_sources', 'test2_sources',
           'test2_sources', 'test2_sources', 'test2_sources']

raynums= ['1000000', '10000000', '100000000',
          '1000000', '10000000', '100000000',
          '1000000', '10000000', '100000000',
          '1000000', '10000000', '100000000',
          '1000000', '10000000', '100000000',
          '1000000', '10000000', '100000000',
          '1000000', '10000000', '100000000',
          '1000000', '10000000', '100000000',
          '1000000', '10000000', '100000000',
          '1000000', '10000000', '100000000']

hydrogen_caseAs= ['False', 'False', 'False',
                  'False', 'False', 'False',
                  'False', 'False', 'False',
                  'False', 'False', 'False',
                  'True', 'True', 'True',
                  'True', 'True', 'True',
                  'True', 'True', 'True',
                  'True', 'True', 'True',
                  'False', 'False', 'False',
                  'False', 'False', 'False']

helium_mfs= ['0.0d0', '0.0d0', '0.0d0',
             '0.0d0', '0.0d0', '0.0d0',
             '0.4d0', '0.4d0', '0.4d0',
             '0.4d0', '0.4d0', '0.4d0',
             '0.0d0', '0.0d0', '0.0d0',
             '0.0d0', '0.0d0', '0.0d0',
             '0.0d0', '0.0d0', '0.0d0',
             '0.0d0', '0.0d0', '0.0d0',
             '0.0d0', '0.0d0', '0.0d0',
             '0.0d0', '0.0d0', '0.0d0']

outdirs= ['../../sphray_output/IT1_N64/r6',           '../../sphray_output/IT1_N64/r7',           '../../sphray_output/IT1_N64/r8',
          '../../sphray_output/IT1_N128/r6',          '../../sphray_output/IT1_N128/r7',          '../../sphray_output/IT1_N128/r8',
          '../../sphray_output/IT1_He_N64/r6',        '../../sphray_output/IT1_He_N64/r7',        '../../sphray_output/IT1_He_N64/r8',
          '../../sphray_output/IT1_He_N128/r6',       '../../sphray_output/IT1_He_N128/r7',       '../../sphray_output/IT1_He_N128/r8',
          '../../sphray_output/IT1_HM01QG_N64/r6',    '../../sphray_output/IT1_HM01QG_N64/r7',    '../../sphray_output/IT1_HM01QG_N64/r8',
          '../../sphray_output/IT1_HM01QG_N128/r6',   '../../sphray_output/IT1_HM01QG_N128/r7',   '../../sphray_output/IT1_HM01QG_N128/r8',
          '../../sphray_output/IT1_HM01QGnd_N64/r6',  '../../sphray_output/IT1_HM01QGnd_N64/r7',  '../../sphray_output/IT1_HM01QGnd_N64/r8',
          '../../sphray_output/IT1_HM01QGnd_N128/r6', '../../sphray_output/IT1_HM01QGnd_N128/r7', '../../sphray_output/IT1_HM01QGnd_N128/r8',
          '../../sphray_output/IT2_N64/r6',           '../../sphray_output/IT2_N64/r7',           '../../sphray_output/IT2_N64/r8',
          '../../sphray_output/IT2_N128/r6',          '../../sphray_output/IT2_N128/r7',          '../../sphray_output/IT2_N128/r8']

outrays= ['10000', '100000', '1000000',
          '10000', '100000', '1000000',
          '10000', '100000', '1000000',
          '10000', '100000', '1000000',
          '10000', '100000', '1000000',
          '10000', '100000', '1000000',
          '10000', '100000', '1000000',
          '10000', '100000', '1000000',
          '10000', '100000', '1000000',
          '10000', '100000', '1000000']

forced_out_files= ['../../rtcp_snapshots/test1_output_times.txt', '../../rtcp_snapshots/test1_output_times.txt', '../../rtcp_snapshots/test1_output_times.txt',
                   '../../rtcp_snapshots/test1_output_times.txt', '../../rtcp_snapshots/test1_output_times.txt', '../../rtcp_snapshots/test1_output_times.txt',
                   '../../rtcp_snapshots/test1_output_times.txt', '../../rtcp_snapshots/test1_output_times.txt', '../../rtcp_snapshots/test1_output_times.txt',
                   '../../rtcp_snapshots/test1_output_times.txt', '../../rtcp_snapshots/test1_output_times.txt', '../../rtcp_snapshots/test1_output_times.txt',
                   '../../rtcp_snapshots/test1_output_times.txt', '../../rtcp_snapshots/test1_output_times.txt', '../../rtcp_snapshots/test1_output_times.txt',
                   '../../rtcp_snapshots/test1_output_times.txt', '../../rtcp_snapshots/test1_output_times.txt', '../../rtcp_snapshots/test1_output_times.txt',
                   '../../rtcp_snapshots/test1_output_times.txt', '../../rtcp_snapshots/test1_output_times.txt', '../../rtcp_snapshots/test1_output_times.txt',
                   '../../rtcp_snapshots/test1_output_times.txt', '../../rtcp_snapshots/test1_output_times.txt', '../../rtcp_snapshots/test1_output_times.txt',
                   '../../rtcp_snapshots/test2_output_times.txt', '../../rtcp_snapshots/test2_output_times.txt', '../../rtcp_snapshots/test2_output_times.txt',
                   '../../rtcp_snapshots/test2_output_times.txt', '../../rtcp_snapshots/test2_output_times.txt', '../../rtcp_snapshots/test2_output_times.txt']



for indx in range(len(raynums)):

    print 'config file: ' + config_files[indx]
    f = open( config_files[indx], 'w' )
    
    for conf, value in zip(config_list, defaults):
        if conf == 'DoTestScenario':
            f.write( conf + ': ' + dotests[indx] + '\n' )                    
        elif conf == 'TestScenario':
            f.write( conf + ': ' + test_names[indx] + '\n' )                    
        elif conf == 'IsoTemp':
            f.write( conf + ': ' + iso_temps[indx] + '\n' )
        elif conf == 'RayDepletion':
            f.write( conf + ': ' + deplete_rays[indx] + '\n' )
        elif conf == 'SpectraFile':
            f.write( conf + ': ' + spectra_files[indx] + '\n' )
        elif conf == 'ParFileBase':
            f.write( conf + ': ' + parbases[indx] + '\n' )
        elif conf == 'SourceFileBase':
            f.write( conf + ': ' + srcbases[indx] + '\n' )
        elif conf == 'ForcedRayNumber':
            f.write( conf + ': ' + raynums[indx] + '\n' )
        elif conf == 'HydrogenCaseA':
            f.write( conf + ': ' + hydrogen_caseAs[indx] + '\n' )
        elif conf == 'He_mf':
            f.write( conf + ': ' + helium_mfs[indx] + '\n' )
        elif conf == 'OutputDir':
            f.write( conf + ': ' + outdirs[indx] + '\n' )
        elif conf == 'IonFracOutRays':
            f.write( conf + ': ' + outrays[indx] + '\n' )
        elif conf == 'ForcedOutFile':
            f.write( conf + ': ' + forced_out_files[indx] + '\n' )
        else:
            f.write( conf + ': ' + value + '\n' )
        
    f.close()
