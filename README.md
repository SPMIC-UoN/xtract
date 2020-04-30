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

---------------------------------------------------------------------

## Citations:


Warrington S, Bryant K, Charquero-Ballester M, Douaud G, Jbabdi S*, Mars R*, Sotiropoulos SN* (Submitted)
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
     xtract -bpx <bedpostX_dir> -out <outputDir> -species <SPECIES> [options]
     xtract -list

     Compulsory arguments:

        -bpx <folder>                          Path to bedpostx folder
        -out <folder>                          Path to output folder
        -species <SPECIES>                     One of HUMAN or MACAQUE

     Optional arguments:
        -list                                  List the tract names used in XTRACT
        -str <file>                            Structures file (format: <tractName> per line OR format: <tractName> [samples=1], 1 means 1000, '#' to skip lines)
        -p <folder>                            Protocols folder (all masks in same standard space) (Default=$FSLDIR/etc/xtract_data/<SPECIES>)
        -stdwarp <std2diff> <diff2std>         Standard2diff and Diff2standard transforms (Default=bedpostx_dir/xfms/{standard2diff,diff2standard})
        -gpu                                   Use GPU version
        -res <mm>                              Output resolution (Default=same as in protocol folders unless '-native' used)
        -ptx_options <options.txt>	           Pass extra probtrackx2 options as a text file to override defaults, e.g. --steplength=0.2 --distthresh=10)

        And EITHER:
        -native                                Run tractography in native (diffusion) space

        OR:
        -ref <refimage> <diff2ref> <ref2diff>  Reference image for running tractography in reference space, Diff2Reference and Reference2Diff transforms

```
---------------------------------------------------------------------

## Running XTRACT:
  XTRACT automatically detects if $SGE_ROOT is set and if so uses FSL_SUB.
  For optimal performance, use the GPU version!!!!

---------------------------------------------------------------------

## Atlases:

- For HUMAN, XTRACT uses the MNI152 standard space in $FSLDIR/etc/standard

- For MACAQUE, XTRACT uses the F99 atlas in Caret - see http://brainvis.wustl.edu/wiki/index.php/Caret:Atlases

  We also provide a copy of the F99 atlas in $FSLDIR/etc/xtract_data/standard/F99. This includes a helper script for registering your own diffusion/structural data to the F99 altas

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

For an example, see $FSLDIR/etc/xtract_data/Human/structureList

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

All the masks above should be in standard space (e.g. MNI152 or F99) if you want to run the same tracking for a collection of subjects.

Next, make a structure file using the format <tractName> <nsamples> per line and call XTRACT using -species <SPECIES> -str <file> -p <folder>, pointing to your new protocols folder 'mytrack'.

---------------------------------------------------------------------

## Visualising results with FSLEYES

The output of XTRACT is a folder that contains tracts in separate folders. We provide a convenient script that can load these tracts (or a subset of the tracts) into FSLEYES using different colours for the different tracts but matching the left/right colours

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
