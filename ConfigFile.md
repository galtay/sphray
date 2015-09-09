# Options #

| **Keyword** | **Type** | **Possible Values / Description** |
|:------------|:---------|:----------------------------------|
| `Verbosity` | Integer  | Level of screen output `{0,1,2,3`}. |
| `DoTestScenario` | Logical  | True, if doing one of the Iliev Tests. |
| `TestScenario` | Character | { `iliev_test1, iliev_test2, iliev_test3, iliev_test4` } |
| `JustInit`  | Logical  | True, stops after initialization. |
| `Comoving`  | Logical  | True, scales quantities using a and h. |
| `IsoTemp`   | Real     | If > 0, fixes temperature of all particles at `IsoTemp`.  `FixSnapTemp` must be `False`. |
| `FixSnapTemp` | Logical  | True, fixes temperature of particles at snapshot values.  Must set `IsoTemp` < 0. |
| `EOStemp`   | Real     | If > 0, fixes temperature of star forming particles at `EOStemp`. |
| `InitxHI`   | Real     | If > 0, initializes all neutral fractions to `InitxHI`. |
| `RayDepletion` | Logical  | True, removes absorbed photons from ray. |
| `IntSeed`   | Integer  | Seed for random number generator. |
| `StaticFieldSimTime` | Real     | For single snapshot runs, determines time to ray trace. |
| `StaticSimTimeUnit` | Character | `{codetime, myr`}                 |
| `InputType` | Integer  | `{1=Gadget Public,2=Gadget CosmoBH,3=Gadget HDF5 OWLS/GIMIC,4=Gadget VBromm, 5=Gadget Public HDF5`} |
| `SnapPath`  | Character | Path to directory containing particle snapshots. |
| `SourcePath` | Character | Path to directory containing source snapshots. |
| `SpectraFile` | Character | Path to spectra file.             |
| `b2cdFile`  | Character | Path to impact parameter -> column depth table file. |
| `AtomicRatesFile` | Character | Path to atomic rates file.        |
| `ParFileBase` | Character | Input particle files = `SnapPath`/`ParFileBase``_``SnapNum``(.FileNumber)``(.hdf5)` |
| `SourceFileBase` | Character | Input source files = `SourcePath`/`SourceFileBase``_``SnapNum` |
| `StartSnapNum` | Integer  | Initial particle snapshot number. |
| `EndSnapNum` | Integer  | Final particle snapshot number.  Set equal to `StartSnapNum` for single snapshot runs. |
| `ParFilesPerSnap` | Integer  | Number of files in each input particle snapshot. |
| `SourceFilesPerSnap` | Integer  | Number of files in each source snapshot.  Must be = 1 for now. |
| `RayScheme` | Character | `{raynum,header`}.  If `raynum`, must be single snapshot run, `ForcedRayNumber` rays are traced.  If `header` the number of rays indicated in the source header files are traced for each snapshot. |
| `RayStats`  | Logical  | True, massive output on each ray traced. |
| `BndryCond` | Integer  | Ray boundary conditions `{-1=Reflective, 0=Transmissive, 1=Periodic`}. |
| `RecRayTol` | Real     | Fraction of particle that needs to recombine to trigger a recombination ray. |
| `RayPhotonTol` | Real     | Fraction of initial photons left in ray to stop updating particles. |
| `WallSampling` | Integer  | Only used if a background source is included `{1=Random, 2=Sobol3D, 3=Sobol2D`}. |
| `OnTheSpotH` | Logical  | True, use 'on the spot' approximation for Hydrogen |
| `OnTheSpotHe` | Logical  | True, use 'on the spot' approximation for Helium |
| `HydrogenCaseA` | Logical  | True, force case A Hydrogen recombination rates even if `OnTheSpotH` = True. |
| `HeliumCaseA` | Logical  | True, force case A Helium recombination rates even if `OnTheSpotHe` = True. |
| `IonTempSolver` | Integer  | `{1=implicit euler, 2=backwards difference`} |
| `Tfloor`    | Real     | Temperature floor.                |
| `Tceiling`  | Real     | Temperature ceiling.              |
| `xfloor`    | Real     | Ionization state floor, between zero and one, should be zero. |
| `xceiling`  | Real     | Ionization state ceiling, between zero and one, should be one. |
| `NeBackground` | Real     | Constant number density of electrons due to metals. |
| `NraysUpdateNoHits` | Integer  | If > 0, forces a zero photoionization rate update on particles that haven't been hit every `NraysUpdateNoHits`. |
| `RecRaysPerSrcRay` | Integer  | Approximate number of recombination rays to trace for each source ray. |
| ` H_mf `    | Real     | Global Hydrogen mass fraction if not specified in particle snapshots. |
| `He_mf`     | Real     | Global Helium mass fraction if not specified in particle snapshots. |
| `OutputDir` | Character | Path to output directory          |
| `OutputFileBase` | Character | Output particle files = `OutputDir`/`OutputFileBase``_``OutputNumber``(.FileNumber)``(.hdf5)` |
| `OutputType` | Integer  | `{1=Binary Gadget with extra fields, 2=HDF5`} |
| `OutputTiming` | Character | `{standard,forced`}.  If `standard`, outputs evenly spaced in time.  If `forced`, outputs specified in `ForcedOutFile`. |
| `NumStdOuts` | Integer  | Number of outputs if `OutputTiming=standard`. |
| `DoInitialOutput` | Logical  | True, Forces an output before any rays are traced. |
| `IonFracOutRays` | Integer  | Produce screen report every `IonFracOutRays` traced. |
| `ForcedOutFile` | Character | Path to file indicating when forced outputs should be made. |
| `ForcedUnits` | Character | `{codetime, myr, mwionfrac`}.  Units for forced output times, code time, Myr, or mass weighted ionization fraction. |
| `PartPerCell` | Integer  | Maximum number of particles per oct-tree leaf |