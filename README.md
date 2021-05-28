## XTRACT - a command-line tool for automated tractography

XTRACT can be used to automatically extract a set of carefully dissected tracts in humans and macaques (other
species to come). It can also be used to define one's own tractography protocols where all the user needs to do is to
define a set of masks in standard space (e.g. MNI152)

The script was written by Saad Jbabdi, Stamatios Sotiropoulos & Shaun Warrington
(based on the autoPtx tool by Marius de Groot - see https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/AutoPtx)

The tractography protocols were created by:

Rogier Mars & Stamatios Sotiropoulos

with help from:
Saad Jbabdi, Kathryn Bryant, Shaun Warrington, Marina Charquero-Ballester, Gwenaelle Douaud

The XTRACT viewer helper script was written by Shaun Warrington

The XTRACT stats helper script was written by Shaun Warrington

---------------------------------------------------------------------

## Citations:


Warrington S, Bryant K, Khrapitchev A, Sallet J, Charquero-Ballester M, Douaud G, Jbabdi S*, Mars R*,
Sotiropoulos SN* (2020) XTRACT - Standardised protocols for automated tractography in the human and
macaque brain. NeuroImage. DOI: 10.1016/j.neuroimage.2020.116923

de Groot M; Vernooij MW. Klein S, Ikram MA, Vos FM, Smith SM, Niessen WJ, Andersson JLR (2013)
Improving alignment in Tract-based spatial statistics: Evaluation and optimization of image registration.
NeuroImage, 76(1), 400-411. DOI: 10.1016/j.neuroimage.2013.03.015


---------------------------------------------------------------------

## Usage:
```
 __  _______ ____      _    ____ _____
 \ \/ /_   _|  _ \    / \  / ___|_   _|
  \  /  | | | |_) |  / _ \| |     | |  
  /  \  | | |  _ <  / ___ \ |___  | |  
 /_/\_\ |_| |_| \_\/_/   \_\____| |_|  


 Usage:
     xtract -bpx <bedpostX_dir> -out <outputDir> -species <SPECIES> [options]
     xtract -bpx <bedpostX_dir> -out <outputDir> -species CUSTOM -str <file> -p <folder> -stdref <reference> [options]
     xtract -list

     Compulsory arguments:

        -bpx <folder>                          Path to bedpostx folder
        -out <folder>                          Path to output folder
        -species <SPECIES>                     One of HUMAN or MACAQUE or CUSTOM

     If -species CUSTOM:
       -str <file>                            Structures file (format: format: <tractName> [samples=1], 1 means 1000, '#' to skip lines)
       -p <folder>                            Protocols folder (all masks in same standard space)
       -stdref <reference>                    Standard space reference image

     Optional arguments:
        -list                                  List the tract names used in XTRACT
        -str <file>                            Structures file (format: <tractName> per line OR format: <tractName> [samples=1], 1 means 1000, '#' to skip lines)
        -p <folder>                            Protocols folder (all masks in same standard space) (Default=$FSLDIR/data/xtract_data/<SPECIES>)
        -stdwarp <std2diff> <diff2std>         Standard2diff and Diff2standard transforms (Default=bedpostx_dir/xfms/{standard2diff,diff2standard})
        -stdref <reference>                    Standard space reference image (Default = $FSLDIR/data/standard/MNI152_T1_1mm [HUMAN], $datadir/standard/F99/mri/struct_brain [MACAQUE])
        -gpu                                   Use GPU version
        -res <mm>                              Output resolution (Default=same as in protocol folders unless '-native' used)
        -ptx_options <options.txt>	            Pass extra probtrackx2 options as a text file to override defaults, e.g. --steplength=0.2 --distthresh=10)

        And EITHER:
        -native                                Run tractography in native (diffusion) space

        OR:
        -ref <refimage> <diff2ref> <ref2diff>  Reference image for running tractography in reference space, Diff2Reference and Reference2Diff transforms

```
---------------------------------------------------------------------

## Running XTRACT:

XTRACT automatically detects if $SGE_ROOT is set and if so uses FSL_SUB. For optimal performance, use the GPU version!!!!

**Outputs of XTRACT**

Under outputDir:

- "commands.txt" - XTRACT processing commands
- "logs" - directory containing the probtrackx log files
- "tracts" - directory continaing tractography results
  - "tractName" - directory per tract, each continaing:
  - "waytotal" - txt file continaing the number of valid streamlines
  - "density.nii.gz" - nifti file containing the fibre probability distribution
  - "density_lenths.nii.gz" - nifti file containing the fibre lengths, i.e. each voxel is the
  average streamline length - this is the "-ompl" probtrackx option
  - "densityNorm.nii.gz" - nifti file continaing the waytotal normalised fibre probability distribution
  (the "density.nii.gz" divided by the total number of valid streamlines)
  - If the protocol calls for reverse-seeding:
    - "tractsInv" - directory continaing the above for the seed-target reversed run
    - "sum_waytotal" and "sum_density.nii.gz" - the summed waytotal and fibre probability distribution
  - If the "-native" option is being used:
    - "masks" - directory
    - "tractName" - directory per tract continaing the native space protocol masks

The primary output is the "densityNorm.nii.gz" file.

**Pre-processing**

Prior to running XTRACT, you must complete the FDT processing pipeline:

1. Brain extraction using BET
2. Susceptibility distortion correction using topup (only if spin-echo fieldmaps have been
  acquired - if you don't have these, skip to step 3)
3. Eddy current distortion and motion correction using eddy
4. Fit the crossing fibre model using bedpostx
5. Registration to standard space (MNI152), see the FDT pipeline

Your data should now be ready to run XTRACT!


---------------------------------------------------------------------

## Atlases:

- For HUMAN, XTRACT uses the MNI152 standard space in $FSLDIR/data/standard

- For MACAQUE, XTRACT uses the F99 atlas in Caret - see http://brainvis.wustl.edu/wiki/index.php/Caret:Atlases

  We also provide a copy of the F99 atlas in $FSLDIR/data/xtract_data/standard/F99. This
  includes a helper script for registering your own diffusion/structural data to the F99 altas

When running XTRACT with the '-species' option, a predefined list of tracts is automatically extracted. Currently the following tracts are available:

| **Tract**   | **Abbreviation** | **XTRACT tractName** |
| -------- | ------------ | ------------ |
| Arcuate Fasciculus | AF | af_l   af_r |
| Acoustic Radiation | AR | ar_l   ar_r |
| Anterior Thalamic Radiation | ATR | atr_l   atr_r |
| Cingulum subsection : Dorsal | CBD | cbd_l   cbd_r |
| Cingulum subsection : Peri-genual | CBP | cbp_l   cbp_r |
| Cingulum subsection : Temporal | CBT | cbt_l   cbt_r |
| Corticospinal Tract | CST | cst_l   cst_r |
| Frontal Aslant | FA | fa_l   fa_r |
| Forceps Major | FMA | fma |
| Forceps Minor | FMI | fmi |
| Fornix | FX | fx_l   fx_r |
| Inferior Longitudinal Fasciculus | ILF | ilf_l   ilf_r |
| Inferior Fronto-Occipital Fasciculus | IFO | ifo_l   ifo_r |
| Middle Cerebellar Peduncle | MCP | mcp |
| Middle Longitudinal Fasciculus | MdLF | mdlf_l   mdlf_r |
| Optic Radiation | OR | or_l or_r |
| Superior Thalamic Radiation | STR | str_l   str_r |
| Superior Longitudinal Fasciculus 1 | SLF1 | slf1_l   slf1_r |
| Superior Longitudinal Fasciculus 2 | SLF2 | slf2_l   slf2_r |
| Superior Longitudinal Fasciculus 3 | SLF3 | slf3_l   slf3_r |
| Anterior Commissure | AC | ac |
| Uncinate Fasciculus | UF | uf_l   uf_r |
| Vertical Occipital Fasciculus | VOF | vof_l   vof_r |

You can run a subset of these tracts by providing a structure text file using the format:

tractName, per line (default number of seeds taken from default structure file)

OR

tractName nsamples, per line

For an example, see $FSLDIR/data/xtract_data/Human/structureList

---------------------------------------------------------------------

## Adding your own tracts:

Suppose you want to create an automated protocol for a tract called 'mytrack'.  

First you need to create a folder called 'mytrack' which you can add e.g. in the protocols folder.

Then create the following NIFTI files (with this exact naming) and copy them into 'mytrack':

**Compulsory**:
- seed.nii.gz : a seed mask

**Optional**:
- stop.nii.gz    : a stop mask if required
- exclude.nii.gz : an exclusion mask if required
- ONE of the following:
  - target.nii.gz  :  a single target mask  
  - target1.nii.gz, target2.nii.gz, etc. : a number of targets, in which case streamlines will be kept if they cross ALL of them
- invert (empty file to indicate that a seed->target and target->seed run will be added and combined)
  if such an option is required a single "target.nii.gz" file is also expected

All the masks above should be in standard space (e.g. MNI152 or F99) if you want to run
the same tracking for a collection of subjects.

Next, make a structure file using the format `<tractName> <nsamples>` per line and call XTRACT
using `-species <SPECIES> -str <file> -p <folder>`, pointing to your new protocols folder 'mytrack'.

---------------------------------------------------------------------

## Running XTRACT with other species:

After developing your own protocols for a given species, you may use the `-species CUSTOM` argument to run XTRACT on any species.
To do so you must provide a standard space brain in addition to specifying the protocol folder and structure file.
e.g.

`xtract -bpx <bedpostX_dir> -out <outputDir> -species CUSTOM -str <file> -p <folder> -stdref <reference> [options]`

---------------------------------------------------------------------

## Visualising results with FSLEYES

The output of XTRACT is a folder that contains tracts in separate folders. We provide a
convenient script that can load these tracts (or a subset of the tracts) into FSLEYES using
different colours for the different tracts but matching the left/right colours

```
 __  _______ ____      _    ____ _____         _                        
 \ \/ /_   _|  _ \    / \  / ___|_   _| __   _(_) _____      _____ _ __
  \  /  | | | |_) |  / _ \| |     | |   \ \ / / |/ _ \ \ /\ / / _ \ '__|
  /  \  | | |  _ <  / ___ \ |___  | |    \ V /| |  __/\ V  V /  __/ |   
 /_/\_\ |_| |_| \_\/_/   \_\____| |_|     \_/ |_|\___| \_/\_/ \___|_|                                                                           


 Usage:
     xtract_viewer -dir <xtractDir> -species HUMAN [options]
     xtract_viewer -dir <xtractDir> -species MACAQUE [options]
     xtract_viewer -dir <xtractDir> -brain <PATH> [options]

     Compulsory arguments:

        -dir <FOLDER>                     Path to XTRACT output folder

        And EITHER:
        -species <SPECIES>                One of HUMAN or MACAQUE

        OR:
        -brain <PATH>                     The brain image to use for the background overlay - must be in the same space as tracts.
                                          Default is the FSL_HCP065_FA map for HUMAN and F99 T1 brain for MACAQUE

     Optional arguments:

        -str STRUCTURE,STRUCTURE,...      Structures (comma separated (default = display all that is found in input folder)

        -thr NUMBER NUMBER                The lower and upper thresholds applied to the tracts for viewing
                                          Default = 0.001 0.1

```

---------------------------------------------------------------------

## Extracting tract-wise summary statistics

A common usage of the XTRACT output is to summarise tracts in terms of simple summary
statistics, such as their volume and microstructural properties (e.g. mean FA). We provide XTRACT
stats to get such summary statistics in a quick and simple way.

You can use XTRACT stats with any modelled diffusion data, e.g. DTI, bedpostx, DKI.

Simply provide; the directory (and basename of files, if any) leading to the diffusion
data of interest, the directory containing the XTRACT output, the warp field (or use 'native'
if tracts are already in diffusion space). If tracts are not in diffusion space, you must also
provide a reference image in diffusion space (e.g. FA map).

e.g. call: `xtract_stats -d /home/DTI/dti_ -xtract /home/xtract -w /home/warp/standard2diff -r /home/DTI/dti_FA`

The output (a .csv file) by default contains the tract volume (mm3) and the mean, median and
standard deviation of the probability, length, FA and MD for each tract.

```
__  _______ ____      _    ____ _____    _        _
\ \/ /_   _|  _ \    / \  / ___|_   _|__| |_ __ _| |_ ___
 \  /  | | | |_) |  / _ \| |     | |/ __| __/ _  | __/ __|
 /  \  | | |  _ <  / ___ \ |___  | |\__ \ || (_| | |_\__ \
/_/\_\ |_| |_| \_\/_/   \_\____| |_||___/\__\__ _|\__|___/


Usage:
    xtract_stats -d <dir_basename> -xtract <XTRACT_dir> -w <xtract2diff> [options]

    Compulsory arguments:

       -d <folder_basename>                   Path to microstructure folder and basename of data (e.g. /home/DTI/dti_)
       -xtract <folder>                       Path to XTRACT output folder
       -w <xtract2diff>                       EITHER XTRACT results to diffusion space transform OR 'native' if tracts are already in diffusion space

    Optional arguments:
       -r <reference>                         If not 'native', provide reference image in diffusion space (e.g. /home/DTI/dti_FA)
       -out <path>                            Output filepath (Default <XTRACT_dir>/stats.csv)
       -str <file>                            Structures file (as in XTRACT) (Default is all tracts under <XTRACT_dir>)
       -thr <float>                           Threshold applied to tract probability map (default = 0.001 = 0.1%)

       -meas <list>                           Comma separated list of features to extract (Default = vol,prob,length,FA,MD - assumes DTI folder has been provided)
                                              vol = tract volume, prob = tract probability, length = tract length
                                              Additional metrics must follow file naming conventions. e.g. for dti_L1 use 'L1'

       -keepfiles                             Keep temporary files

```

---------------------------------------------------------------------

## Generating Connectivity Blueprints with xtract_blueprint

xtract_blueprint calculates the "connectivity blueprint" that represents how different grey matter areas are connected to a set of white matter bundles (Mars et al. 2018 eLife, Warrington et al. 2020 NeuroImage). Formally, this is a seeds x tracts matrix where each row describes how each (cortical) seed location is connected to major white matter fibre bundles, and each column describes the cortical terminations of each of the white matter fibre bundles. xtract_blueprint can, in principle, build a connectivity blueprint for any brain (e.g. different species), at any resolution (i.e. number of surface vertices) and for any region (whole whole or a given ROI, a lobe for example). It inherently supports surface input files (for instance a GIFTI representing the WM/GM boundary interface).

The script was written by Shaun Warrington (University of Nottingham), Saad Jbabdi (Oxford University) and Stamatios Sotiropoulos (University of Nottingham).

---------------------------------------------------------------------

**Citations:**

Mars R, Sotiropoulos SN, Passingham RE, Sallet J, Verhagen L, Khrapitchev AA, Sibson N, Jbabdi S (2018) Whole brain comparative anatomy using connectivity blueprints. eLife. DOI: 10.7554/eLife.35237

Warrington S, Bryant K, Khrapitchev A, Sallet J, Charquero-Ballester M, Douaud G, Jbabdi S*, Mars R*, Sotiropoulos SN* (2020) XTRACT - Standardised protocols for automated tractography in the human and macaque brain. NeuroImage. DOI: 10.1016/j.neuroimage.2020.116923

---------------------------------------------------------------------
**Usage:**
```

__  _______ ____      _    ____ _____ _     _                       _       _
\ \/ /_   _|  _ \    / \  / ___|_   _| |__ | |_   _  ___ _ __  _ __(_)_ __ | |_
 \  /  | | | |_) |  / _ \| |     | | | '_ \| | | | |/ _ \ '_ \| '__| | '_ \| __|
 /  \  | | |  _ <  / ___ \ |___  | | | |_) | | |_| |  __/ |_) | |  | | | | | |_
/_/\_\ |_| |_| \_\/_/   \_\____| |_| |_.__/|_|\__,_|\___| .__/|_|  |_|_| |_|\__|
                                                        |_|


Usage:
    xtract_blueprint -bpx <folder> -xtract <folder> -seeds <list> -warps <ref> <xtract2diff> <diff2xtract> -out <folder> [options]

    Compulsory arguments:

       -bpx       <folder>                          Path to bedpostx folder
       -out       <folder>                          Path to output folder
       -xtract    <folder>                          Path to xtract folder
       -seeds     <list>                            Comma separated list of seeds for which a blueprint is requested (e.g. left and right cortex in standard space)
       -warps     <ref> <xtract2diff> <diff2xtract> Standard space reference image and transforms between xtract space and diffusion space

    Optional arguments:
       -stage                                       What to run. 1:matrix2, 2:blueprint, all:everythng (default)
       -gpu                                         Use GPU version
       -savetxt                                     Save blueprint as txt file (nseed by ntracts) instead of CIFTI
       -prefix    <string>                          Specify a prefix for the final blueprint filename (e.g. <prefix>_BP.LR.dscalar.nii)

       -rois  <list>                                Comma separated list (1 per seed): ROIs (gifti) to restrict seeding (e.g. medial wall masks)
       -stops     <stop.txt>                        Text file containing line separated list
       -wtstops   <wtstop.txt>                      Text file containing line separated list
       -tract_list                                  Comma separated list of tracts to include (default = all found under -xtract <folder>)

       -thr                                         Threshold applied to XTRACT tracts prior to blueprint calculation (default = 0.001, i.e. 0.1% probability).
       -nsamples                                    Number of samples per seed used in tractography (default = 1000)
       -res       <mm>                              Resolution of matrix2 output (Default = 3 mm)
       -ptx_options <options.txt>                   Pass extra probtrackx2 options as a text file to override defaults

   Example (recommended) usage:
      xtract_blueprint -bpx sub001/dMRI.bedpostx -out sub001/blueprint -xtract sub001/xtract -seeds sub001/l.white.surf.gii,sub001/r.white.surf.gii \
           -rois sub001/l.medwall.shape.gii,sub001/r.medwall.shape.gii -warps MNI152_brain.nii.gz sub001/xtract2diff.nii.gz sub001/diff2xtract.nii.gz -gpu \


```

---------------------------------------------------------------------

## Running xtract_blueprint:

In order to use xtract_blueprint, you need to have run xtract first. To construct the connectivity blueprints, xtract_blueprint expects the same warp fields as used in the running of xtract.

xtract_blueprint currently supports the construction of connectivity blueprints using surface (GIFTI) files only. As such, you must provide GIFTI surface files containing the white-grey matter boundary surface data. (We plan to extend this to support volume seeding in the future.)

**Required input:**
-	bedpostx folder - crossing fibre modelled diffusion data (expects to find nodif_brain_mask)
-	xtract folder   - xtract tract folder (the output from xtract)
-	seed 			      - the comma separated seed masks to use in tractography (e.g. the white-grey matter boundary surfaces), e.g. `L.white.surf.gii,R.white.surf.gii`
-	warps			      - a reference standard space image and warps to and from the native diffusion and standard spaces, e.g. `MNI152.nii.gz standard2diff.nii.gz diff2standard.nii.gz`

Note: if running whole-brain (recommended), you must provide the left seed first, as exampled.
xtract_blueprint uses the seed GIFTI header information to build the CIFTI blueprint. Ensure that your structure is set in the seed GIFTI: either CortexLeft or CortexRight.

**Running modes:**
- Stage 1 only		- only run seed-based tractography
- Stage 2 only		- only run blueprint processing (requires xtract output and tractography from stage 1)
- All			        - runs both stage 1 and stage 2 processing

xtract_blueprint is capable of GPU acceleration (`'-gpu'` flag) and detects $SGE_ROOT to work with fsl_sub. If using the CPU version, expect tractography to take many hours.


**Tractography details and options:**

Tractography is performed for each seed separately. The resultant connectivity matrices (fdt_matrix) are stacked in order to construct a whole-brain connectivity blueprint. You may also provide a single hemisphere if you wish.

Optionally, you may also provide stop (stop tracking at locations given by this mask file) and wtstop (allow propagation within mask but terminate on exit) masks. Stop is typically the pial surface. wtstop is typically subcortical structures and/or the white surface. These should be specified as line separated text files. e.g. `-seeds <l.white.surf.gii,r.white.surf.gii> -stop stop.txt -wtstop wtstop.txt`

Spatial resolution: by default tractography will be ran and stored using a resolution of 3 mm. This may be adjusted using the `'-res'` argument. Note: if required, xtract_blueprint will resample the xtract tracts. Warning: connectivity matrices are very large and require a lot of memory to handle - 3 mm is usually sufficient for the adult human brain.

Additional probtrackx options may also be supplied. Simply add the probtrackx arguments to a text file and direct xtract_blueprint to this using the `'-ptx_options'` argument.

Connectivity blueprints are primarily concerned with the connectivity of the cortex to white matter tracts. As such, we offer the option the mask out the medial wall. To do so, provide a single medial wall mask per supplied seed: e.g. `-seeds l.white.surf.gii,r.white.surf.gii -rois l.roi,r.roi`. By default, the medial wall is included in the calculation of the connectivity blueprint: we recommend the use of the medial wall mask to prevent this.

The '-roi' argument may be used to restrict the blueprint to any region of interest, not just to exclude the medial wall. For example, you may provide an ROI restricting the blueprint to the temporal or frontal lobe.

If you wish to use a stop/wtstop surface mask, you must ensure that the number of vertices matches the seed mask. This means that, if you are providing a seed mask and medial wall mask to xtract_blueprint, and wish to provide a surface stop mask, you must convert the stop mask to asc, restricting the data points to the medial wall mask, e.g.:

	${FSLDIR}/bin/surf2surf -i stop.L.surf.gii -o stop.L.asc --outputtype=ASCII --values=l.roi.shape.gii
	${FSLDIR}/bin/surf2surf -i stop.R.surf.gii -o stop.R.asc --outputtype=ASCII --values=r.roi.shape.gii

Then supply "stop.L.asc" and "stop.R.asc" in a text file to xtract_blueprint using the `'-stop'` argument. This conversion is automatically performed for the seed mask in xtract_blueprint if a medial wall mask is supplied.


**Which tracts are included?**

Connectivity blueprints may be constructed using the provided XTRACT tracts or using your own. By default, xtract_blueprint will use all tracts it finds under the xtract folder. You can specify a subset, or you own tracts, by providing a comma separated list of tracts using the `'-tracts <str,str,str>'` argument.

Certain tracts, e.g. the Middle Cerebellar Peduncle (MCP), do not project to the cortex. As such, they should be disregarded when interpreting the connectivity blueprint, or excluded from its construction.


**Only interested in the connectivity of a specific area?**

In many cases, the connectivity to a particular lobe, e.g. temporal or frontal, is of interest. You can use xtract_blueprint to obtain a connectivity blueprint for such a region:

1. Define a binary mask which contains the region of interest as a shape.gii or func.gii file
2. Select the tracts of interest: in all likelihood, only a subset of XTRACT's tracts will project to the ROI
3. Supply the whole white matter surface file along with the ROI to xtract_blueprint, e.g. for the temporal lobe

  xtract_blueprint -bpx sub001/dMRI.bedpostx -out sub001/blueprint -xtract sub001/xtract \
  -warps MNI152_brain.nii.gz sub001/xtract2diff.nii.gz sub001/diff2xtract.nii.gz -gpu \
  -seeds sub001/l.white.surf.gii -rois sub001/l.temporal_lobe.shape.gii -tract_list af_l,ilf_l,ifo_l,mdlf_l


**Outputs of xtract_blueprint**

xtract_blueprint will create an output directory specified by the `'-out'` argument. This will contain any log and command files along with a sub-directory per seed. Each sub-directory contains the resultant connectivity matrix from stage 1. The connectivity blueprint will be saved in the parent output directory (a CIFTI dscalar.nii file).

Under outputDir:
- Stage 1 (matrix2 tractography) output
    - omatrix2:
        - "ptx_commands.txt" - the probtrackx commands for tractography
        - "`<seed>`" - sub-directory containing tractography results for each seed supplied, each containing:
            - "coords_fdt_matrix2" - the coordinates of the seed mask
            - "lookup_tractspace_fdt_matrix.nii.gz" - the target lookup space
            - "tract_space_coords_for_fdt_matrix2" - the coordinates of the target mask
            - "probtrackx.log" - the probtrackx log file
            - "fdt_paths.nii.gz" - a NIFTI file containing the generated streamlines
            - "fdt_matrix.dot" - the sparse format connectivity matrix (used to calculated the blueprint)
            - "waytotal" - txt file continaing the number of valid streamlines
- Stage 2 (blueprint calculation) output
    - "bp_commands.txt" - the blueprint calculation commands
    - "BP.`<L,R,LR>`.dscalar.nii" - CIFTI file containing the whole-brain connectivity blueprint - if running both hemispheres
- "logs" - sub-directory containing job scheduler log files for both stages

Alternatively, the `'-savetxt'` option may be used to override this. In this case, two txt files will be saved: the first (BP.<L,R,LR>.txt) will be an n_seed by n_tracts array containing the blueprint; the second (tract_order.`<L,R,LR>`.txt) is an n_tracts by 1 array containing the tract order in which the blueprint is structured. Note: no CIFTI file will be generated.
