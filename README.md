# openFoam-Overset-SpeedUp


This repository provides a method (small code patches, parameters, and decomposition methodology) to achieve significant **speed and precision improvements** with respect to the current implementation of the overset (Chimera) code of the ESI-Version of OpenFOAM. This is achieved by modifying the search function, in particular the search boxes and their parallelization and by pushing the interpolation layer of the background mesh further away from the walls of the component.

The reported tests were conducted on OpenFOAM-v2012 using overPimpleFyMFoam, but the method should work almost *as is* for latter releases and for other solvers.

A quick and dirty presentation on the obtained speedups and impact on precision when compared to AMI is given in the [accompanying presentation](https://github.com/louisgag/openFoam-Overset-SpeedUp/raw/main/oversetComparisonGagnon.pdf).

## How to (assuming you start from a working overset case)

### 1. compile both the modified decomposePar app and overset libs
Load openFoam-v2012, go into the inverseDistancePushFront folder and run wmake
Load openFoam-v2106, go to the myDecomposePar folder and run wmake (may run on v2012 as well but I haven't tested because I use a different machine to decompose and to run the case)

### 2. decompose what will never be in the interpolation zone separately from what will be
Use the myDecomposePar code in order to split the domain into two regions and first decompose them independently to later join them, thus ensuring that the processors that are out of your region of interest never disturb the search process by, for example, creating huge search boxes.
*You need to calculate an appropriate number of processors for each zone (i.e. each overset component) that will lead to a balanced number of cells per processors.*
Make a copy of your case and proceed as such:
#### a) On the copied case (where you have myDecomposePar and probably v2106):
`topoSet -dict system/topoSetDict_c0_in # set zones that will correspond to the separately decomposed zones; ensure to delete all existing zones before doing this (you don't have to delete them on the original case)`

`splitMeshRegions -cellZonesFileOnly cellZones -overwrite`

`cd system;for i in "c0" "c0_in" "bl0" "bl1" "bl2";do echo $i;cd $i;ln -sT ../decomposeParDict ./decomposeParDict;cd ..;done;cd .. # this is only an example where there is one background zone which gets divided into two regions (where interpolation can occur, c0_in, and where it can't occur, c0, along with the component meshes, bl1,2,3`

`myDecomposePar -allRegions -copyZero -cellDist -force > log.decomp 2>&1  ## you need to defined a regionProperties file in the constant folder where you list the names of your regions and say they are all fluid regions, -cellDist is required otherwise globalProcIds will not be written... copyZero allows you to save the time required to write the fields, which you don't need yet`
#### b) on the original case (where you will run the case and don't need myDecomposePar):
copy constant/globalProcIds that was created by myDecomposePar and run

`decomposePar -force -cellDist -decomposeParDict system/decomposeParDictMan  ## dry-run will not work, -cellDist allows you to see if it works (as a hint go into paraview make a copy of your search boxes using a paraview source and ensure no cell from your outter regions even touch the border of your search boxes.`

The decomposeParDictMan should be defined as such:
```
numberOfSubdomains 1024; # your total number of processors

method manual;
manualCoeffs
{
    dataFile    "globalProcIds";
}
```

### 3. include the overset lib in your controlDict file:
  `libs            (overset oversetIDpushfront);`

### 4. choose the method and an appropriate voxel size in you fvSchemes file:
  ```
oversetInterpolation
{
      method inverseDistancePushFront;
      searchBox           (-1 -1 -0.55)(1 1 0.55); // this box should be **entirely** included in the inner decomposed domain
      voxelSize 0.008; // the voxel size size determines how fine the search for donor and acceptor cells will be, too fine will slow down your simulation and too large will make it instable. Take into account that too coarse voxels may also lead to donor cells being outside of your inner decomposed box.
      nPushFront 14; // push interpolation front away from hole by n layers (0 to disable)
      layerRelax 1;
}
```
### 5. run the case
e.g. `mpirun -n 1024 omplace -nt 1 -c 0-125:st=3,127 overPimpleDyMFoam -parallel -fileHandler collated > log.pimpleFoam.$(date +%s) 2>&1;`
### 6. Tips
- if you work with pimpleFoam, use a *final* relaxation of 0.9 on all fields and equations, it will make the pressure oscillations almost completly disappear and have practically no impact on the accuracy of your results (see the plots of slides 5 and 6 on the accompanying pdf).
  - why? partly because when pimple decides it converged, it still does one last iteration and therefore relaxing there is not terrible.
- if you run into division by 0 (floating point error) during the calculation of the wall distance, decrease the tolerance on yPsi (1e-10 for example).
- unless your mesh is huge use single-double precision: `export WM_PRECISION_OPTION=SPDP;`
- both `leastSquaresPushFront` and `trackingInverseDistancePushFront` that come along with this patch were not thoroughly tested and did not give satisfactory results, use **inverseDistancePushFront** unless you know what you are doing or have a lot of time to play around.
- this code is released to help others trying to get a faster overset search, but is also not particularly clean or well commented.
- the myDeomposePar was originally intented to [make AMI faster](http://louisgagnon.com/articlesScientifiques/LouisGagnonOFW2020SplashTalk_vFinal.pdf), by allowing to decompose separately each AMI-zone, but it only ended up making everything slower. It thus only became useful when attacking the overset problem, where it really makes things much faster.
- the [voxelSize contribution](https://develop.openfoam.com/Development/openfoam/-/issues/1683) came from [Nicolas Edh](https://develop.openfoam.com/nicolasedh) who graciously accepted that I include it in this code.
