

## XTRACT - a command-line tool for automated tractography

XTRACT can be used to automatically extract a set of carefully dissected tracts in humans and macaques (other 
species to come). It can also be used to define one's own tractography protocols where all the user needs to do is to 
define a set of masks in standard space (e.g. MNI)

The script was written by Saad Jbabdi & Stamatios Sotiropoulos
(based on the autoPtx tool by Marius de Groot - see https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/AutoPtx)

The tractography protocols were created by:

Rogier Mars & Stamatios Sotiropoulos

with help from:
Saad Jbabdi, Kathryn Bryant, Shaun Warrington, Marina Charquero-Ballester, Gwenaelle Douaud


---------------------------------------------------------------------

## Citations:


Warrington S, Bryant K, Charquero-Ballester M, Douaud G, Jbabdi S*, Mars R*, Sotiropoulos SN* (in prep.)
Standardised protocols for automated tractography and connectivity blueprints in the human and macaque brain.

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
    xtract -bpx <bedpostX_dir> -out <outputDir> -str <structuresFile> -p <protocolsFolder> [options]
    xtract -bpx <bedpostX_dir> -out <outputDir> -species HUMAN [options]
    xtract -bpx <bedpostX_dir> -out <outputDir> -species MACAQUE [options]

    Compulsory arguments:

       -bpx <folder>                     Path to bedpostx folder
       -out <folder>                     Path to output folder
       
       And EITHER:
       -str <file>                       Structures file (format: <tractName> [samples=1], 1 means 1000, '#' to skip lines)
       -p   <folder>                     Protocols folder (all masks in same standard space)

       Or:
       -species <SPECIES>                One of HUMAN or MACAQUE

    Optional arguments:

       -stdwarp <std2diff> <diff2std>    Standard2diff and Diff2standard transforms (Default=bedpostx_dir/xfms/{standard2diff,diff2standard}) 
       -gpu                              Use GPU version 
       -native                           Run tractography in native (diffusion) space
       -res <mm>                         Output resolution (Default=same as in protocol folders unless '-native' used)

```
---------------------------------------------------------------------

## Running XTRACT:
  XTRACT automatically detects if $SGE_ROOT is set and if so uses FSL_SUB 
  For optimal performance, use the GPU version!!!! 

---------------------------------------------------------------------

## Atlases:

- For HUMAN, XTRACT uses the MNI152 standard space in $FSLDIR/etc/standard

- For MACAQUE, XTRACT uses the F99 atlas in Caret - see http://brainvis.wustl.edu/wiki/index.php/Caret:Atlases
  
  We also provide a copy of the F99 atlas in $FSLDIR/etc/xtract_data/standard/F99
  This includes a helper script for registering your own diffusion/structural data to the F99 altas

When running XTRACT with the '-species' option, a predefined list of tracts is automatically extracted. Currently the following tracts are available:

| **Tract**   | **Abbreviation** |
| --------| ------------ |
| Arcuate Fasciculus | AF |
| Acoustic Radiation | AR |
| Anterior Thalamic Radiation | ATR |
| Cingulum subsection : Dorsal | CBD |
| Cingulum subsection : Parahippocampal | CBP |
| Cingulum subsection : Temporal | CBT |
| Corticospinal Tract | CST |
| Frontal Aslant | FA |
| Forceps Major | FMA |
| Forceps Minor | FMI |
| Fornix | FX |
| Inferior Longitudinal Fasciculus | ILF |
| Inferior Fronto-Occipital Fasciculus | IFO |
| Middle Cerebellar Peduncle | MCP |
| Middle Longitudinal Fasciculuc | MdLF |
| Optic Radiation | OR |
| Superior Thalamic Radiation | STR |
| Superior Longitudinal Fasciculus 1 | SLF1 |
| Superior Longitudinal Fasciculus 2 | SLF2 |
| Superior Longitudinal Fasciculus 3 | SLF3 |
| Anterior Commissure | AC |
| Uncinate Fasciculus | UF |
| Vertical Occipital Fasciculus | VOF |



---------------------------------------------------------------------

## Adding your own tracts:

Suppose you want to create an automated protocol for a tract called 'mytrack'.  

First you need to create a folder called 'mytrack' which you can add e.g. in the protocols folder. 

Then create the following NIFTI files (with this exact naming) and copy them into mytrack:

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

All the masks above should be in standard space (e.g. MNI152 or F99) if you want to run the same tracking for a collection of subjects.


---------------------------------------------------------------------

## Visualising results with FSLEYES

The output of XTRACT is a folder that contrains tracts in separate folders. We provide a convenient script that can load these tracts (or a subset of the tracts) into FSLEYES using different colours for the different tracts but matching the left/right colours

```
 __  _______ ____      _    ____ _____         _                        
 \ \/ /_   _|  _ \    / \  / ___|_   _| __   _(_) _____      _____ _ __ 
  \  /  | | | |_) |  / _ \| |     | |   \ \ / / |/ _ \ \ /\ / / _ \ '__|
  /  \  | | |  _ <  / ___ \ |___  | |    \ V /| |  __/\ V  V /  __/ |   
 /_/\_\ |_| |_| \_\/_/   \_\____| |_|     \_/ |_|\___| \_/\_/ \___|_|                                                                           
                                                             

Usage:
    xtract_viewer -dir <xtractDir> [options]

    Compulsory arguments:

       -dir FOLDER                       Path to XTRACT output folder

    Optional arguments:

       -str STRUCTURE,STRUCTURE,...     Structures (comma separated (default = display all that is found in input folder)

       -thr NUMBER NUMBER                The lower and upper thresholds applied to the tracts for viewing
                                         Default = 0.001 0.1

       -brain                            The brain image to use for the background overlay - must be in the same space as tracts.
                                         Default is the FSL_HCP065_FA map

```
