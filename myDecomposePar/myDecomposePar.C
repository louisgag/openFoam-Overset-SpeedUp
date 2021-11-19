/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    decomposePar

Group
    grpParallelUtilities

Description
    Automatically decomposes a mesh and fields of a case for parallel
    execution of OpenFOAM.

Usage
    \b decomposePar [OPTIONS]

    Options:
      - \par -allRegions
        Decompose all regions in regionProperties. Does not check for
        existence of processor*.

      - \par -case \<dir\>
        Specify case directory to use (instead of the cwd).

      - \par -cellDist
        Write the cell distribution as a labelList, for use with 'manual'
        decomposition method and as a VTK or volScalarField for visualization.

      - \par -constant
        Include the 'constant/' dir in the times list.

      - \par -copyUniform
        Copy any \a uniform directories too.

      - \par -copyZero
        Copy \a 0 directory to processor* rather than decompose the fields.

      - \par -debug-switch \<name=val\>
        Specify the value of a registered debug switch. Default is 1
        if the value is omitted. (Can be used multiple times)

      - \par -decomposeParDict \<file\>
        Use specified file for decomposePar dictionary.

      - \par -dry-run
        Test without writing the decomposition. Changes -cellDist to
        only write VTK output.

      - \par -fields
        Use existing geometry decomposition and convert fields only.

      - \par fileHandler \<handler\>
        Override the file handler type.

      - \par -force
        Remove any existing \a processor subdirectories before decomposing the
        geometry.

      - \par -ifRequired
        Only decompose the geometry if the number of domains has changed from a
        previous decomposition. No \a processor subdirectories will be removed
        unless the \a -force option is also specified. This option can be used
        to avoid redundant geometry decomposition (eg, in scripts), but should
        be used with caution when the underlying (serial) geometry or the
        decomposition method etc. have been changed between decompositions.

      - \par -info-switch \<name=val\>
        Specify the value of a registered info switch. Default is 1
        if the value is omitted. (Can be used multiple times)

      - \par -latestTime
        Select the latest time.

      - \par -lib \<name\>
        Additional library or library list to load (can be used multiple times).

      - \par -noFunctionObjects
        Do not execute function objects.

      - \par -noSets
        Skip decomposing cellSets, faceSets, pointSets.

      - \par -noZero
        Exclude the \a 0 dir from the times list.

      - \par -opt-switch \<name=val\>
        Specify the value of a registered optimisation switch (int/bool).
        Default is 1 if the value is omitted. (Can be used multiple times)

      - \par -region \<regionName\>
        Decompose named region. Does not check for existence of processor*.

      - \par -time \<ranges\>
        Override controlDict settings and decompose selected times. Does not
        re-decompose the mesh i.e. does not handle moving mesh or changing
        mesh cases. Eg, ':10,20 40:70 1000:', 'none', etc.

      - \par -verbose
        Additional verbosity.

      - \par -doc
        Display documentation in browser.

      - \par -doc-source
        Display source code in browser.

      - \par -help
        Display short help and exit.

      - \par -help-man
        Display full help (manpage format) and exit.

      - \par -help-notes
        Display help notes (description) and exit.

      - \par -help-full
        Display full help and exit.

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "fvCFD.H"
#include "IOobjectList.H"
#include "domainDecomposition.H"
#include "domainDecompositionDryRun.H"
#include "labelIOField.H"
#include "labelFieldIOField.H"
#include "scalarIOField.H"
#include "scalarFieldIOField.H"
#include "vectorIOField.H"
#include "vectorFieldIOField.H"
#include "sphericalTensorIOField.H"
#include "sphericalTensorFieldIOField.H"
#include "symmTensorIOField.H"
#include "symmTensorFieldIOField.H"
#include "tensorIOField.H"
#include "tensorFieldIOField.H"
#include "pointFields.H"
#include "regionProperties.H"

#include "readFields.H"
#include "fvFieldDecomposer.H"
#include "pointFieldDecomposer.H"
#include "lagrangianFieldDecomposer.H"
#include "decompositionModel.H"

#include "faCFD.H"
#include "emptyFaPatch.H"
#include "faMeshDecomposition.H"
#include "faFieldDecomposer.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Read proc addressing at specific instance.
// Uses polyMesh/fvMesh meshSubDir by default
autoPtr<labelIOList> procAddressing
(
    const fvMesh& procMesh,
    const word& name,
    const word& instance,
    const word& local = fvMesh::meshSubDir
)
{
    return autoPtr<labelIOList>::New
    (
        IOobject
        (
            name,
            instance,
            local,
            procMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false  // do not register
        )
    );
}


// Read proc addressing at specific instance.
// Uses the finiteArea meshSubDir
autoPtr<labelIOList> faProcAddressing
(
    const fvMesh& procMesh,
    const word& name,
    const word& instance,
    const word& local = faMesh::meshSubDir
)
{
    return procAddressing(procMesh, name, instance, local);
}


// Return cached or read proc addressing from facesInstance
const labelIOList& procAddressing
(
    const PtrList<fvMesh>& procMeshList,
    const label proci,
    const word& name,
    PtrList<labelIOList>& procAddressingList
)
{
    const fvMesh& procMesh = procMeshList[proci];

    if (!procAddressingList.set(proci))
    {
        procAddressingList.set
        (
            proci,
            procAddressing(procMesh, name, procMesh.facesInstance())
        );
    }
    return procAddressingList[proci];
}


void decomposeUniform
(
    const bool copyUniform,
    const domainDecomposition& mesh,
    const Time& processorDb,
    const word& regionDir = word::null
)
{
    const Time& runTime = mesh.time();

    // Any uniform data to copy/link?
    const fileName uniformDir(regionDir/"uniform");

    if (fileHandler().isDir(runTime.timePath()/uniformDir))
    {
        Info<< "Detected additional non-decomposed files in "
            << runTime.timePath()/uniformDir
            << endl;

        const fileName timePath =
            fileHandler().filePath(processorDb.timePath());

        // If no fields have been decomposed the destination
        // directory will not have been created so make sure.
        mkDir(timePath);

        if (copyUniform || mesh.distributed())
        {
            if (!fileHandler().exists(timePath/uniformDir))
            {
                fileHandler().cp
                (
                    runTime.timePath()/uniformDir,
                    timePath/uniformDir
                );
            }
        }
        else
        {
            // Link with relative paths
            string parentPath = string("..")/"..";

            if (regionDir != word::null)
            {
                parentPath = parentPath/"..";
            }

            fileName currentDir(cwd());
            chDir(timePath);

            if (!fileHandler().exists(uniformDir))
            {
                fileHandler().ln
                (
                    parentPath/runTime.timeName()/uniformDir,
                    uniformDir
                );
            }
            chDir(currentDir);
        }
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Decompose a mesh and fields of a case for parallel execution"
    );

    argList::noParallel();
    argList::addOption
    (
        "decomposeParDict",
        "file",
        "Use specified file for decomposePar dictionary"
    );

    #include "addAllRegionOptions.H"

    argList::addBoolOption
    (
        "dry-run",
        "Test without writing the decomposition. "
        "Changes -cellDist to only write VTK output."
    );
    argList::addBoolOption
    (
        "verbose",
        "Additional verbosity"
    );
    argList::addOption
    (
        "domains",
        "N",
        "Override numberOfSubdomains (-dry-run only)",
        true  // Advanced option
    );
    argList::addOption
    (
        "method",
        "name",
        "Override decomposition method (-dry-run only)",
        true  // Advanced option
    );

    argList::addBoolOption
    (
        "cellDist",
        "Write cell distribution as a labelList - for use with 'manual' "
        "decomposition method and as a volScalarField for visualization."
    );
    argList::addBoolOption
    (
        "copyZero",
        "Copy 0/ directory to processor*/ rather than decompose the fields"
    );
    argList::addBoolOption
    (
        "copyUniform",
        "Copy any uniform/ directories too"
    );
    argList::addBoolOption
    (
        "fields",
        "Use existing geometry decomposition and convert fields only"
    );
    argList::addBoolOption
    (
        "noSets",
        "Skip decomposing cellSets, faceSets, pointSets"
    );
    argList::addBoolOption
    (
        "force",
        "Remove existing processor*/ subdirs before decomposing the geometry"
    );
    argList::addBoolOption
    (
        "ifRequired",
        "Only decompose geometry if the number of domains has changed"
    );

    // Allow explicit -constant, have zero from time range
    timeSelector::addOptions(true, false);  // constant(true), zero(false)

    #include "setRootCase.H"

    const bool dryrun           = args.found("dry-run");
    const bool writeCellDist    = args.found("cellDist");
    const bool verbose          = args.found("verbose");
    const bool mergeCellDist    = 1; // create a global cellDist for a region-based decomposition

    // Most of these are ignored for dry-run (not triggered anywhere)
    const bool copyZero         = args.found("copyZero");
    const bool copyUniform      = args.found("copyUniform");
    const bool decomposeSets    = !args.found("noSets");
    const bool decomposeIfRequired = args.found("ifRequired");

    bool decomposeFieldsOnly = args.found("fields");
    bool forceOverwrite      = args.found("force");


    // Set time from database
    #include "createTime.H"

    // Allow override of time (unless dry-run)
    instantList times;
    if (dryrun)
    {
        Info<< "\ndry-run: ignoring -copy*, -fields, -force, time selection"
            << nl;
    }
    else
    {
        times = timeSelector::selectIfPresent(runTime, args);
    }

    // Allow override of decomposeParDict location
    fileName decompDictFile(args.get<fileName>("decomposeParDict", ""));
    if (!decompDictFile.empty() && !decompDictFile.isAbsolute())
    {
        decompDictFile = runTime.globalPath()/decompDictFile;
    }

    // Get region names
    #include "getAllRegionOptions.H"

    const bool optRegions =
        (regionNames.size() != 1 || regionNames[0] != polyMesh::defaultRegion);

    if (regionNames.size() == 1 && regionNames[0] != polyMesh::defaultRegion)
    {
        Info<< "Using region: " << regionNames[0] << nl << endl;
    }
    
    labelList globalProcIds; // cannot be out of scope for upcoming loop...
    label procOffset = 0; // offset to apply to processors of current region
    if (mergeCellDist)
    {
        Info << "mergeCellDist arg given, creating global list of processors (with offset, todo)" << endl;
    }

    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        const word& regionDir =
        (
            regionName == polyMesh::defaultRegion ? word::null : regionName
        );

        if (dryrun)
        {
            Info<< "dry-run: decomposing mesh " << regionName << nl << nl
                << "Create mesh..." << flush;
                
            Info<< "LG regioni: " << regioni << endl;

            domainDecompositionDryRun decompTest
            (
                IOobject
                (
                    regionName,
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                decompDictFile,
                args.getOrDefault<label>("domains", 0),
                args.getOrDefault<word>("method", word::null)
            );

            decompTest.execute(writeCellDist, verbose);
            continue;
        }

        Info<< "\n\nDecomposing mesh";
        if (!regionDir.empty())
        {
            Info<< ' ' << regionName;
        }
        Info<< nl << endl;

        // Determine the existing processor count directly
        label nProcs = fileHandler().nProcs(runTime.path(), regionDir);

        // Get requested numberOfSubdomains directly from the dictionary.
        // Note: have no mesh yet so cannot use decompositionModel::New
        const label nDomains = decompositionMethod::nDomains
        (
            IOdictionary
            (
                IOobject::selectIO
                (
                    IOobject
                    (
                        decompositionModel::canonicalName,
                        runTime.time().system(),
                        regionDir,  // region (if non-default)
                        runTime,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false       // do not register
                    ),
                    decompDictFile
                )
            )
        );

        // Give file handler a chance to determine the output directory
        const_cast<fileOperation&>(fileHandler()).setNProcs(nDomains);

        if (decomposeFieldsOnly)
        {
            // Sanity check on previously decomposed case
            if (nProcs != nDomains)
            {
                FatalErrorInFunction
                    << "Specified -fields, but the case was decomposed with "
                    << nProcs << " domains"
                    << nl
                    << "instead of " << nDomains
                    << " domains as specified in decomposeParDict" << nl
                    << exit(FatalError);
            }
        }
        else if (nProcs)
        {
            bool procDirsProblem = true;

            if (decomposeIfRequired && nProcs == nDomains)
            {
                // We can reuse the decomposition
                decomposeFieldsOnly = true;
                procDirsProblem = false;
                forceOverwrite = false;

                Info<< "Using existing processor directories" << nl;
            }

            if (optRegions)
            {
                procDirsProblem = false;
                forceOverwrite = false;
            }

            if (forceOverwrite)
            {
                Info<< "Removing " << nProcs
                    << " existing processor directories" << endl;

                // Remove existing processors directory
                fileNameList dirs
                (
                    fileHandler().readDir
                    (
                        runTime.path(),
                        fileName::Type::DIRECTORY
                    )
                );
                forAllReverse(dirs, diri)
                {
                    const fileName& d = dirs[diri];

                    label proci = -1;

                    if
                    (
                        d.starts_with("processor")
                     &&
                        (
                            // Collated is "processors"
                            d[9] == 's'

                            // Uncollated has integer(s) after 'processor'
                         || Foam::read(d.substr(9), proci)
                        )
                    )
                    {
                        if (fileHandler().exists(d))
                        {
                            fileHandler().rmDir(d);
                        }
                    }
                }

                procDirsProblem = false;
            }

            if (procDirsProblem)
            {
                FatalErrorInFunction
                    << "Case is already decomposed with " << nProcs
                    << " domains, use the -force option or manually" << nl
                    << "remove processor directories before decomposing. e.g.,"
                    << nl
                    << "    rm -rf " << runTime.path().c_str() << "/processor*"
                    << nl
                    << exit(FatalError);
            }
        }

        Info<< "Create mesh" << endl;
        domainDecomposition mesh
        (
            IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            decompDictFile
        );

        // Decompose the mesh
        if (!decomposeFieldsOnly)
        {
            mesh.decomposeMesh();

            mesh.writeDecomposition(decomposeSets);

            if (writeCellDist)
            {
                const labelList& procIds = mesh.cellToProc();

                // Write decomposition for visualization
                mesh.writeVolField("cellDist");
                //TBD: mesh.writeVTK("cellDist");

                // Write decomposition as labelList for use with 'manual'
                // decomposition method.
                labelIOList cellDecomposition
                (
                    IOobject
                    (
                        "cellDecomposition",
                        mesh.facesInstance(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    procIds
                );
                cellDecomposition.write();

                Info<< nl << "Wrote decomposition to "
                    << runTime.relativePath(cellDecomposition.objectPath())
                    << " for use in manual decomposition." << endl;
/*
         cellZoneMesh newCellZones
         (
             IOobject
             (
                 zoneFile,
                 mesh.facesInstance(),
                 polyMesh::meshSubDir,
                 mesh,
                 IOobject::MUST_READ,
                 IOobject::NO_WRITE,
                 false
             ),
             mesh
         );
                    */
                    
    //~ labelList regIds;
                    
     labelIOList parentIds
     (
         IOobject
         (
             "cellRegionAddressing",
             mesh.facesInstance(),
             polyMesh::meshSubDir,
             mesh,
             IOobject::MUST_READ,
             IOobject::NO_WRITE,
             false
         )//,
         //~ regIds
     );
     
//     Info << "LG: parentIds.size() = " << parentIds.size() << "; max(parentIds) = " << max(parentIds) << endl;
//     Info << "LG: globalProcIds.size() = " << globalProcIds.size() <<  endl;

     
     if (max(parentIds) >= globalProcIds.size())
     {
	     globalProcIds.resize(max(parentIds)+1);
     }
                 forAll(parentIds, celli)
                    {
//                        Info << "LG celli = " << celli <<  "; parentIds[celli] = " << parentIds[celli] << "; procIds[celli] = " << procIds[celli] << ";" << endl;
                        globalProcIds[parentIds[celli]] = procIds[celli]+procOffset;
                    }   
		 Info << "LG: decomposed region " << regionName << " into " << max(procIds)+1
			 << " processors assigned to ranks from " << procOffset << " to "
			 << procOffset + max(procIds) << endl;
		 procOffset += max(procIds)+1;
            }

            fileHandler().flush();
        }


        if (copyZero)
        {
            // Copy the 0 directory into each of the processor directories
            fileName prevTimePath;
            for (label proci = 0; proci < mesh.nProcs(); ++proci)
            {
                Time processorDb
                (
                    Time::controlDictName,
                    args.rootPath(),
                    args.caseName()/("processor" + Foam::name(proci))
                );
                processorDb.setTime(runTime);

                if (fileHandler().isDir(runTime.timePath()))
                {
                    // Get corresponding directory name (to handle processors/)
                    const fileName timePath
                    (
                        fileHandler().objectPath
                        (
                            IOobject
                            (
                                "",
                                processorDb.timeName(),
                                processorDb
                            ),
                            word::null
                        )
                    );

                    if (timePath != prevTimePath)
                    {
                        Info<< "Processor " << proci
                            << ": copying " << runTime.timePath() << nl
                            << " to " << timePath << endl;
                        fileHandler().cp(runTime.timePath(), timePath);

                        prevTimePath = timePath;
                    }
                }
            }
        }
        else
        {
            // Decompose the field files

            // Cached processor meshes and maps. These are only preserved if
            // running with multiple times.
            PtrList<Time> processorDbList(mesh.nProcs());
            PtrList<fvMesh> procMeshList(mesh.nProcs());
            PtrList<labelIOList> faceProcAddressingList(mesh.nProcs());
            PtrList<labelIOList> cellProcAddressingList(mesh.nProcs());
            PtrList<labelIOList> boundaryProcAddressingList(mesh.nProcs());
            PtrList<fvFieldDecomposer> fieldDecomposerList(mesh.nProcs());
            PtrList<labelIOList> pointProcAddressingList(mesh.nProcs());
            PtrList<pointFieldDecomposer> pointFieldDecomposerList
            (
                mesh.nProcs()
            );


            // Loop over all times
            forAll(times, timeI)
            {
                runTime.setTime(times[timeI], timeI);

                Info<< "Time = " << runTime.timeName() << endl;

                // Search for list of objects for this time
                IOobjectList objects(mesh, runTime.timeName());


                // Construct the vol fields
                // ~~~~~~~~~~~~~~~~~~~~~~~~
                PtrList<volScalarField> volScalarFields;
                readFields(mesh, objects, volScalarFields, false);
                PtrList<volVectorField> volVectorFields;
                readFields(mesh, objects, volVectorFields, false);
                PtrList<volSphericalTensorField> volSphericalTensorFields;
                readFields(mesh, objects, volSphericalTensorFields, false);
                PtrList<volSymmTensorField> volSymmTensorFields;
                readFields(mesh, objects, volSymmTensorFields, false);
                PtrList<volTensorField> volTensorFields;
                readFields(mesh, objects, volTensorFields, false);


                // Construct the dimensioned fields
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                PtrList<DimensionedField<scalar, volMesh>> dimScalarFields;
                readFields(mesh, objects, dimScalarFields);
                PtrList<DimensionedField<vector, volMesh>> dimVectorFields;
                readFields(mesh, objects, dimVectorFields);
                PtrList<DimensionedField<sphericalTensor, volMesh>>
                    dimSphericalTensorFields;
                readFields(mesh, objects, dimSphericalTensorFields);
                PtrList<DimensionedField<symmTensor, volMesh>>
                    dimSymmTensorFields;
                readFields(mesh, objects, dimSymmTensorFields);
                PtrList<DimensionedField<tensor, volMesh>> dimTensorFields;
                readFields(mesh, objects, dimTensorFields);


                // Construct the surface fields
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                PtrList<surfaceScalarField> surfaceScalarFields;
                readFields(mesh, objects, surfaceScalarFields, false);
                PtrList<surfaceVectorField> surfaceVectorFields;
                readFields(mesh, objects, surfaceVectorFields, false);
                PtrList<surfaceSphericalTensorField>
                    surfaceSphericalTensorFields;
                readFields(mesh, objects, surfaceSphericalTensorFields, false);
                PtrList<surfaceSymmTensorField> surfaceSymmTensorFields;
                readFields(mesh, objects, surfaceSymmTensorFields, false);
                PtrList<surfaceTensorField> surfaceTensorFields;
                readFields(mesh, objects, surfaceTensorFields, false);


                // Construct the point fields
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~
                const pointMesh& pMesh = pointMesh::New(mesh);

                PtrList<pointScalarField> pointScalarFields;
                readFields(pMesh, objects, pointScalarFields, false);
                PtrList<pointVectorField> pointVectorFields;
                readFields(pMesh, objects, pointVectorFields, false);
                PtrList<pointSphericalTensorField> pointSphericalTensorFields;
                readFields(pMesh, objects, pointSphericalTensorFields, false);
                PtrList<pointSymmTensorField> pointSymmTensorFields;
                readFields(pMesh, objects, pointSymmTensorFields, false);
                PtrList<pointTensorField> pointTensorFields;
                readFields(pMesh, objects, pointTensorFields, false);


                // Construct the Lagrangian fields
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                fileNameList cloudDirs
                (
                    fileHandler().readDir
                    (
                        runTime.timePath()/cloud::prefix,
                        fileName::DIRECTORY
                    )
                );

                // Particles
                PtrList<Cloud<indexedParticle>> lagrangianPositions
                (
                    cloudDirs.size()
                );
                // Particles per cell
                PtrList<List<SLList<indexedParticle*>*>> cellParticles
                (
                    cloudDirs.size()
                );

                PtrList<PtrList<labelIOField>> lagrangianLabelFields
                (
                    cloudDirs.size()
                );
                PtrList<PtrList<labelFieldCompactIOField>>
                lagrangianLabelFieldFields
                (
                    cloudDirs.size()
                );
                PtrList<PtrList<scalarIOField>> lagrangianScalarFields
                (
                    cloudDirs.size()
                );
                PtrList<PtrList<scalarFieldCompactIOField>>
                lagrangianScalarFieldFields
                (
                    cloudDirs.size()
                );
                PtrList<PtrList<vectorIOField>> lagrangianVectorFields
                (
                    cloudDirs.size()
                );
                PtrList<PtrList<vectorFieldCompactIOField>>
                lagrangianVectorFieldFields
                (
                    cloudDirs.size()
                );
                PtrList<PtrList<sphericalTensorIOField>>
                lagrangianSphericalTensorFields
                (
                    cloudDirs.size()
                );
                PtrList<PtrList<sphericalTensorFieldCompactIOField>>
                    lagrangianSphericalTensorFieldFields(cloudDirs.size());
                PtrList<PtrList<symmTensorIOField>> lagrangianSymmTensorFields
                (
                    cloudDirs.size()
                );
                PtrList<PtrList<symmTensorFieldCompactIOField>>
                lagrangianSymmTensorFieldFields
                (
                    cloudDirs.size()
                );
                PtrList<PtrList<tensorIOField>> lagrangianTensorFields
                (
                    cloudDirs.size()
                );
                PtrList<PtrList<tensorFieldCompactIOField>>
                lagrangianTensorFieldFields
                (
                    cloudDirs.size()
                );

                label cloudI = 0;

                for (const fileName& cloudDir : cloudDirs)
                {
                    IOobjectList cloudObjects
                    (
                        mesh,
                        runTime.timeName(),
                        cloud::prefix/cloudDir,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    );

                    // Note: look up "positions" for backwards compatibility
                    if
                    (
                        cloudObjects.found("coordinates")
                     || cloudObjects.found("positions")
                    )
                    {
                        // Read lagrangian particles
                        // ~~~~~~~~~~~~~~~~~~~~~~~~~

                        Info<< "Identified lagrangian data set: "
                            << cloudDir << endl;

                        lagrangianPositions.set
                        (
                            cloudI,
                            new Cloud<indexedParticle>
                            (
                                mesh,
                                cloudDir,
                                false
                            )
                        );


                        // Sort particles per cell
                        // ~~~~~~~~~~~~~~~~~~~~~~~

                        cellParticles.set
                        (
                            cloudI,
                            new List<SLList<indexedParticle*>*>
                            (
                                mesh.nCells(),
                                static_cast<SLList<indexedParticle*>*>(nullptr)
                            )
                        );

                        label i = 0;

                        for (indexedParticle& p : lagrangianPositions[cloudI])
                        {
                            p.index() = i++;

                            label celli = p.cell();

                            // Check
                            if (celli < 0 || celli >= mesh.nCells())
                            {
                                FatalErrorInFunction
                                    << "Illegal cell number " << celli
                                    << " for particle with index "
                                    << p.index()
                                    << " at position "
                                    << p.position() << nl
                                    << "Cell number should be between 0 and "
                                    << mesh.nCells()-1 << nl
                                    << "On this mesh the particle should"
                                    << " be in cell "
                                    << mesh.findCell(p.position())
                                    << exit(FatalError);
                            }

                            if (!cellParticles[cloudI][celli])
                            {
                                cellParticles[cloudI][celli] =
                                    new SLList<indexedParticle*>();
                            }

                            cellParticles[cloudI][celli]->append(&p);
                        }

                        // Read fields
                        // ~~~~~~~~~~~

                        IOobjectList lagrangianObjects
                        (
                            mesh,
                            runTime.timeName(),
                            cloud::prefix/cloudDirs[cloudI],
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE,
                            false
                        );

                        lagrangianFieldDecomposer::readFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianLabelFields
                        );

                        lagrangianFieldDecomposer::readFieldFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianLabelFieldFields
                        );

                        lagrangianFieldDecomposer::readFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianScalarFields
                        );

                        lagrangianFieldDecomposer::readFieldFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianScalarFieldFields
                        );

                        lagrangianFieldDecomposer::readFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianVectorFields
                        );

                        lagrangianFieldDecomposer::readFieldFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianVectorFieldFields
                        );

                        lagrangianFieldDecomposer::readFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianSphericalTensorFields
                        );

                        lagrangianFieldDecomposer::readFieldFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianSphericalTensorFieldFields
                        );

                        lagrangianFieldDecomposer::readFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianSymmTensorFields
                        );

                        lagrangianFieldDecomposer::readFieldFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianSymmTensorFieldFields
                        );

                        lagrangianFieldDecomposer::readFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianTensorFields
                        );

                        lagrangianFieldDecomposer::readFieldFields
                        (
                            cloudI,
                            lagrangianObjects,
                            lagrangianTensorFieldFields
                        );

                        cloudI++;
                    }
                }

                lagrangianPositions.setSize(cloudI);
                cellParticles.setSize(cloudI);
                lagrangianLabelFields.setSize(cloudI);
                lagrangianLabelFieldFields.setSize(cloudI);
                lagrangianScalarFields.setSize(cloudI);
                lagrangianScalarFieldFields.setSize(cloudI);
                lagrangianVectorFields.setSize(cloudI);
                lagrangianVectorFieldFields.setSize(cloudI);
                lagrangianSphericalTensorFields.setSize(cloudI);
                lagrangianSphericalTensorFieldFields.setSize(cloudI);
                lagrangianSymmTensorFields.setSize(cloudI);
                lagrangianSymmTensorFieldFields.setSize(cloudI);
                lagrangianTensorFields.setSize(cloudI);
                lagrangianTensorFieldFields.setSize(cloudI);

                Info<< endl;

                // split the fields over processors
                for (label proci = 0; proci < mesh.nProcs(); ++proci)
                {
                    Info<< "Processor " << proci << ": field transfer" << endl;

                    // open the database
                    if (!processorDbList.set(proci))
                    {
                        processorDbList.set
                        (
                            proci,
                            new Time
                            (
                                Time::controlDictName,
                                args.rootPath(),
                                args.caseName()
                              / ("processor" + Foam::name(proci))
                            )
                        );
                    }
                    Time& processorDb = processorDbList[proci];


                    processorDb.setTime(runTime);

                    // read the mesh
                    if (!procMeshList.set(proci))
                    {
                        procMeshList.set
                        (
                            proci,
                            new fvMesh
                            (
                                IOobject
                                (
                                    regionName,
                                    processorDb.timeName(),
                                    processorDb
                                )
                            )
                        );
                    }
                    const fvMesh& procMesh = procMeshList[proci];

                    const labelIOList& faceProcAddressing = procAddressing
                    (
                        procMeshList,
                        proci,
                        "faceProcAddressing",
                        faceProcAddressingList
                    );

                    const labelIOList& cellProcAddressing = procAddressing
                    (
                        procMeshList,
                        proci,
                        "cellProcAddressing",
                        cellProcAddressingList
                    );

                    const labelIOList& boundaryProcAddressing = procAddressing
                    (
                        procMeshList,
                        proci,
                        "boundaryProcAddressing",
                        boundaryProcAddressingList
                    );


                    // FV fields: volume, surface, internal
                    {
                        if (!fieldDecomposerList.set(proci))
                        {
                            fieldDecomposerList.set
                            (
                                proci,
                                new fvFieldDecomposer
                                (
                                    mesh,
                                    procMesh,
                                    faceProcAddressing,
                                    cellProcAddressing,
                                    boundaryProcAddressing
                                )
                            );
                        }
                        const fvFieldDecomposer& fieldDecomposer =
                            fieldDecomposerList[proci];

                        // vol fields
                        fieldDecomposer.decomposeFields(volScalarFields);
                        fieldDecomposer.decomposeFields(volVectorFields);
                        fieldDecomposer.decomposeFields
                        (
                            volSphericalTensorFields
                        );
                        fieldDecomposer.decomposeFields(volSymmTensorFields);
                        fieldDecomposer.decomposeFields(volTensorFields);

                        // surface fields
                        fieldDecomposer.decomposeFields(surfaceScalarFields);
                        fieldDecomposer.decomposeFields(surfaceVectorFields);
                        fieldDecomposer.decomposeFields
                        (
                            surfaceSphericalTensorFields
                        );
                        fieldDecomposer.decomposeFields
                        (
                            surfaceSymmTensorFields
                        );
                        fieldDecomposer.decomposeFields(surfaceTensorFields);

                        // internal fields
                        fieldDecomposer.decomposeFields(dimScalarFields);
                        fieldDecomposer.decomposeFields(dimVectorFields);
                        fieldDecomposer.decomposeFields(dimSphericalTensorFields);
                        fieldDecomposer.decomposeFields(dimSymmTensorFields);
                        fieldDecomposer.decomposeFields(dimTensorFields);

                        if (times.size() == 1)
                        {
                            // Clear cached decomposer
                            fieldDecomposerList.set(proci, nullptr);
                        }
                    }


                    // Point fields
                    if
                    (
                        pointScalarFields.size()
                     || pointVectorFields.size()
                     || pointSphericalTensorFields.size()
                     || pointSymmTensorFields.size()
                     || pointTensorFields.size()
                    )
                    {
                        const labelIOList& pointProcAddressing = procAddressing
                        (
                            procMeshList,
                            proci,
                            "pointProcAddressing",
                            pointProcAddressingList
                        );

                        const pointMesh& procPMesh = pointMesh::New(procMesh);

                        if (!pointFieldDecomposerList.set(proci))
                        {
                            pointFieldDecomposerList.set
                            (
                                proci,
                                new pointFieldDecomposer
                                (
                                    pMesh,
                                    procPMesh,
                                    pointProcAddressing,
                                    boundaryProcAddressing
                                )
                            );
                        }
                        const pointFieldDecomposer& pointDecomposer =
                            pointFieldDecomposerList[proci];

                        pointDecomposer.decomposeFields(pointScalarFields);
                        pointDecomposer.decomposeFields(pointVectorFields);
                        pointDecomposer.decomposeFields
                        (
                            pointSphericalTensorFields
                        );
                        pointDecomposer.decomposeFields(pointSymmTensorFields);
                        pointDecomposer.decomposeFields(pointTensorFields);


                        if (times.size() == 1)
                        {
                            pointProcAddressingList.set(proci, nullptr);
                            pointFieldDecomposerList.set(proci, nullptr);
                        }
                    }


                    // If there is lagrangian data write it out
                    forAll(lagrangianPositions, cloudI)
                    {
                        if (lagrangianPositions[cloudI].size())
                        {
                            lagrangianFieldDecomposer fieldDecomposer
                            (
                                mesh,
                                procMesh,
                                faceProcAddressing,
                                cellProcAddressing,
                                cloudDirs[cloudI],
                                lagrangianPositions[cloudI],
                                cellParticles[cloudI]
                            );

                            // Lagrangian fields
                            {
                                fieldDecomposer.decomposeFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianLabelFields[cloudI]
                                );
                                fieldDecomposer.decomposeFieldFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianLabelFieldFields[cloudI]
                                );
                                fieldDecomposer.decomposeFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianScalarFields[cloudI]
                                );
                                fieldDecomposer.decomposeFieldFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianScalarFieldFields[cloudI]
                                );
                                fieldDecomposer.decomposeFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianVectorFields[cloudI]
                                );
                                fieldDecomposer.decomposeFieldFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianVectorFieldFields[cloudI]
                                );
                                fieldDecomposer.decomposeFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianSphericalTensorFields[cloudI]
                                );
                                fieldDecomposer.decomposeFieldFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianSphericalTensorFieldFields[cloudI]
                                );
                                fieldDecomposer.decomposeFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianSymmTensorFields[cloudI]
                                );
                                fieldDecomposer.decomposeFieldFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianSymmTensorFieldFields[cloudI]
                                );
                                fieldDecomposer.decomposeFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianTensorFields[cloudI]
                                );
                                fieldDecomposer.decomposeFieldFields
                                (
                                    cloudDirs[cloudI],
                                    lagrangianTensorFieldFields[cloudI]
                                );
                            }
                        }
                    }

                    // Decompose the "uniform" directory in the time region
                    // directory
                    decomposeUniform(copyUniform, mesh, processorDb, regionDir);

                    // For a multi-region case, also decompose the "uniform"
                    // directory in the time directory
                    if (regionNames.size() > 1 && regioni == 0)
                    {
                        decomposeUniform(copyUniform, mesh, processorDb);
                    }

                    // We have cached all the constant mesh data for the current
                    // processor. This is only important if running with
                    // multiple times, otherwise it is just extra storage.
                    if (times.size() == 1)
                    {
                        boundaryProcAddressingList.set(proci, nullptr);
                        cellProcAddressingList.set(proci, nullptr);
                        faceProcAddressingList.set(proci, nullptr);
                        procMeshList.set(proci, nullptr);
                        processorDbList.set(proci, nullptr);
                    }
                }

                // Finite area mesh and field decomposition

                IOobject faMeshBoundaryIOobj
                (
                    "faBoundary",
                    mesh.time().findInstance
                    (
                        mesh.dbDir()/polyMesh::meshSubDir,
                        "boundary"
                    ),
                    faMesh::meshSubDir,
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false  // not registered
                );


                if (faMeshBoundaryIOobj.typeHeaderOk<faBoundaryMesh>(true))
                {
                    Info<< "\nFinite area mesh decomposition" << endl;

                    // Always based on the volume decomposition!
                    faMeshDecomposition aMesh
                    (
                        mesh,
                        mesh.nProcs(),
                        mesh.model()
                    );

                    aMesh.decomposeMesh();
                    aMesh.writeDecomposition();


                    // Construct the area fields
                    // ~~~~~~~~~~~~~~~~~~~~~~~~
                    PtrList<areaScalarField> areaScalarFields;
                    readFields(aMesh, objects, areaScalarFields);

                    PtrList<areaVectorField> areaVectorFields;
                    readFields(aMesh, objects, areaVectorFields);

                    PtrList<areaSphericalTensorField> areaSphericalTensorFields;
                    readFields(aMesh, objects, areaSphericalTensorFields);

                    PtrList<areaSymmTensorField> areaSymmTensorFields;
                    readFields(aMesh, objects, areaSymmTensorFields);

                    PtrList<areaTensorField> areaTensorFields;
                    readFields(aMesh, objects, areaTensorFields);


                    // Construct the edge fields
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    PtrList<edgeScalarField> edgeScalarFields;
                    readFields(aMesh, objects, edgeScalarFields);

                    const label nAreaFields =
                    (
                        areaScalarFields.size()
                      + areaVectorFields.size()
                      + areaSphericalTensorFields.size()
                      + areaSymmTensorFields.size()
                      + areaTensorFields.size()
                      + edgeScalarFields.size()
                    );

                    Info<< endl;
                    Info<< "Finite area field transfer: "
                        << nAreaFields << " fields" << endl;

                    // Split the fields over processors
                    for
                    (
                        label proci = 0;
                        nAreaFields && proci < mesh.nProcs();
                        ++proci
                    )
                    {
                        Info<< "    Processor " << proci << endl;

                        // open the database
                        Time processorDb
                        (
                            Time::controlDictName,
                            args.rootPath(),
                            args.caseName()
                          / ("processor" + Foam::name(proci))
                        );

                        processorDb.setTime(runTime);

                        // Read the mesh
                        fvMesh procFvMesh
                        (
                            IOobject
                            (
                                regionName,
                                processorDb.timeName(),
                                processorDb
                            )
                        );

                        faMesh procMesh(procFvMesh);

                        // // Does not work.  HJ, 15/Aug/2017
                        // const labelIOList& faceProcAddressing =
                        //     procAddressing
                        //     (
                        //         procMeshList,
                        //         proci,
                        //         "faceProcAddressing",
                        //         faceProcAddressingList
                        //     );

                        // const labelIOList& boundaryProcAddressing =
                        //     procAddressing
                        //     (
                        //         procMeshList,
                        //         proci,
                        //         "boundaryProcAddressing",
                        //         boundaryProcAddressingList
                        //     );

                        // Addressing from faMesh (not polyMesh) meshSubDir

                        autoPtr<labelIOList> tfaceProcAddr =
                            faProcAddressing
                            (
                                procFvMesh,
                                "faceProcAddressing",
                                runTime.constant()
                            );
                        auto& faceProcAddressing = *tfaceProcAddr;

                        autoPtr<labelIOList> tboundaryProcAddr =
                            faProcAddressing
                            (
                                procFvMesh,
                                "boundaryProcAddressing",
                                runTime.constant()
                            );
                        auto& boundaryProcAddressing = *tboundaryProcAddr;

                        autoPtr<labelIOList> tedgeProcAddr =
                            faProcAddressing
                            (
                                procFvMesh,
                                "edgeProcAddressing",
                                runTime.constant()
                            );
                        const auto& edgeProcAddressing = *tedgeProcAddr;

                        faFieldDecomposer fieldDecomposer
                        (
                            aMesh,
                            procMesh,
                            edgeProcAddressing,
                            faceProcAddressing,
                            boundaryProcAddressing
                        );

                        fieldDecomposer.decomposeFields(areaScalarFields);
                        fieldDecomposer.decomposeFields(areaVectorFields);
                        fieldDecomposer.decomposeFields
                        (
                            areaSphericalTensorFields
                        );
                        fieldDecomposer.decomposeFields(areaSymmTensorFields);
                        fieldDecomposer.decomposeFields(areaTensorFields);

                        fieldDecomposer.decomposeFields(edgeScalarFields);
                    }
                }
            }
        }
    }

            if (writeCellDist)
            {
    Info<< "LG: globalProcIds[0] = " << globalProcIds[0] << endl;
    Info<< "LG: globalProcIds.size() = " << globalProcIds.size() << endl;
    Info<< "LG: globalProcIds[globalProcIds.size()] = " << globalProcIds[globalProcIds.size()] << endl;
    Info<< "LG: globalProcIds[globalProcIds.size()+1] = " << globalProcIds[globalProcIds.size()+1] << endl;
/*
   forAll (globalProcIds,celli)
    {
	   Info << "LG globalProcIds[celli] = " << globalProcIds[celli] << "; celli = " << celli << endl;
    }
*/

                // Write decomposition as labelList for use with 'manual'
                // decomposition method.
                labelIOList globalCellDecomposition
                (
                    IOobject
                    (
                        "globalProcIds",
                        runTime.constant(),
                        runTime,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                   globalProcIds
                );
                globalCellDecomposition.write();

                Info<< nl << "Wrote global decomposition to "
                    << globalCellDecomposition.objectPath()
                    << " for use in manual decomposition." << endl;
	    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
