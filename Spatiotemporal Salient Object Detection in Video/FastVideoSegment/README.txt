Fast video segmentation v1.2
================================================================================

A. Papazoglou and Vittorio Ferrari
a.papazoglou@sms.ed.ac.uk

This software was developed under 64-bit Linux with Matlab R2012b/R2013a.There
is no guarantee it will run on other operating systems or Matlab versions
(though it probably will).

If you use this software for academic research, please cite [1].

All software with the exception of the contents of the folder 'Code/External/'
are released as free software for non-commercial applications under the terms of
the GNU General Public License version 3 (GPLv3), see 
<http://www.gnu.org/licenses/>, which is also noted in each individual file. For
commercial applications, contact the author for licensing options.

This software uses components from [2,3,4,5,6], which are included in
'Code/External'. Please check their individual release licences.

The Matlab wrapper functions for the aforementioned components are released
under the BSD license, as also noted in each individual file. Specifically, the
wrapper files published under the BSD license are:
    - Code/External/maxflow/maxflow_mex_optimisedWrapper.cpp
    - Code/External/SLIC/SLIC_mex.cpp
    - Code/External/sundaramECCV10_ldof_GPU_mex.cpp
    - compile.m
    
If you find any bugs or have any comments, please e-mail A. Papazoglou.



Introduction
--------------------------------------------------------------------------------
This release implements the fast video object segmentation technique presented
in [1]. The algorithm requires optical flow between subsequent frames, as well
per-frame superpixel over-segmentation as input. To quickly test our method we
have included the optical flow estimation method of [2,3], and the superpixel
over-segmentation methods [4,5]. If you wish to use alternative methods instead,
refer to the “Using alternative optical flow/superpixels methods” section.

In this release archive you will find:

    - The source code implementation of [1].
    - The optical flow methods of [2,3], as provided by their authors.
    - The superpixel methods [4,5] as provided by their authors.
    - An implementation of the maxflow algorithm [6] as provided by the authors.
    - Matlab mex wrappers for the code of [3,4,5,6].
    - Compiled mex files for 64-bit Linux.
    - Two short sequences from the YouTube-Objects dataset [7].



Quick start
--------------------------------------------------------------------------------
Let us assume that <DIR> is the directory where you uncompressed the release
archive. In the Matlab prompt, type:

    cd <DIR>
    compile   #compiles all the nescessary mex files
    demo

Note: The release archive already contains precompiled mex files for 64-bit
Linux, so you only need to run the compile.m script if you are on another
platform (e.g. 32-bit Linux or Windows).

The demo.m script segments two test shots from YouTube-Objects, located in
‘<DIR>/Data/inputs/animals/’. In detail, demo.m does the following:

    1) Computes the optical flow between subsequent frames. The default method
       used is [2].
    2) Computes per-frame superpixel over-segmentations. The default method used
       is [5].
    3) Computes the video object segmentation as described in [1].

The outputs (optical flow, superpixels, segmentations, visuals) will be stored
in ‘<DIR>/Data/outputs/animals/’. '<DIR>/Data/model outputs/animals/' contains
precomputed segmentations (using the default parameters). If you are using the
software correctly, your outputs will match the model outputs.

For further details, check the comments in demo.m directly. If you intend to use
this software as part of a larger pipeline, or to test it on your own dataset
you should read the next sections as well.


Note: If - and only if - you want to use the GPU accelerated optical flow method
of [3], then you should use 64-bit Linux, as it is only provided in binary form
by the authors. This requires CUDA 5.0 to be properly installed on your system,
and it should also be on the default library path.



options struct
--------------------------------------------------------------------------------
Most wrapper functions require an options struct that specifies where the input
data are stored, preprocessing methods to be used etc. If you intend to use the
code beyond the scope of the demo, you need to create your own options struct.
The necessary fields that the algorithm expects are:

options.infolder
    The full path to the folder where the input frames are stored. The code
    expects the framesto be stored in .jpg format, and be named sequentially
    using 8 digits (e.g 00000001.jpg, 00000002.jpg etc). If you want to use
    another file format or naming convention you need to modify
    ‘<DIR>/Code/Loaders/readFrame.m’. Note that the folder can
    contain multiple videos/shots, as long as all the frames follow the naming
    convention.

options.outfolder
    The full path to the folder where the outputs will be stored. If it does not
    already exist, the code will create it automatically.

options.ranges
    A matlab array of length S+1, containing the number of the first frame of
    each shot (assuming that there are S shots are present in options.infolder).
    The last element of the array should be equal to the number of frames + 1.

options.visualise
    This field should be either true or false, and denotes whether the code
    should create videos of the various processing stages (for visual inspection
    of the results and/or debugging).

options.vocal
    This field should be either true or false, and denotes whether the code
    should print status and timing messages during execution.

Some optional fields used by the preprocessing functions:

options.flowmethod
    This field denotes the optical flow method to be used by
    computeOpticalFlow.m and loadFlow.m. The valid arguments are ‘broxPAMI2011’
    and ‘sundaramECCV2010’. ‘sundaramECCV2010’ is typically much faster, but
    requires CUDA 5.0 to be installed as well as a CUDA capable GPU.

options.superpixels
    This field denotes the superpixel oversegmentation method to be used by
    computeSuperpixels.m and loadSuperpixels.m. The valid arguments are
    ‘Turbopixels’ and ‘SLIC’. Turbopixels is the method used for the results in
    [1], while we use SLIC in the demo.m here.
Both methods work equally well and ‘SLIC’ is much faster.



params struct
--------------------------------------------------------------------------------
Besides the options struct, the ‘videoRapidSegment.m’ function expects a params
struct containing numerical values for weights etc. For details, we refer to
[1]. The required fields are:

params.maxIterations
    The maximum number of iterations for the refinement step.
    The value used in [1] was 4, as we typically do not see much change beyond
    the 3rd or 4th iteration.

params.fadeout
    The temporal weighting factor for appearance learning ([1], eq. 10).

params.foregroundMixtures, params.backgroundMixtures
    The number of mixture components for foreground/background
    appearance modelling.

params.spatialWeight, params.temporalWeight
    The weights of the spatial and temporal smoothness potentials([1], eq. 6).

The function getDefaultParams.m provides some sane default values for these
parameters (i.e. they seem to give reasonable segmentations for a large number
of videos). These are the parameters used to produce the large-scale
YouTube-objects experiments in [1].  Note, however, that these default values
have not been optimised on any dataset. If you have ground-truth annotations
for the dataset that you intend to apply the algorithm to, we suggest that you
optimise these values on that dataset.



Using alternative optical flow/superpixels methods
--------------------------------------------------------------------------------
The code is optical flow/superpixel method agnostic, so you could use your own
methods if you wish. The only requirement is that you use the same data format:

Optical flow:
    The optical flow should be stored in a cell array (N-1)x1, where N is the
    number of frames in the shot. Each cell should contain a HxWx2 int16 matrix
    that contains the displacement values, where H is the height of the frame
    and W is the width. flow{ frame }( :, :, 1 ) should contain the displacement
    over the height axis, while  flow{ frame }( :, :, 2 ) should contain
    the displacement over the width axis.

Superpixels:
    The superpixels should be stored in a cell array Nx1. Each cell should
    contain a HxW uint16 matrix that depicts the superpixel label/ID for each
    pixel in the frame. The superpixel labels need to be > 0.
    
Both of these cell need to be passed as part of the input struct to the
videoRapidSegment.m function, as data.flow and data.superpixels respectively.
For an example, see demo.m.


Running the algorithm on the SegTrack dataset
--------------------------------------------------------------------------------
SegTrack videos are drastically different to the YouTube-Objects videos, in size
-both frame size and typical object size- as well as quality. If you want to
evaluate the performance of the algorithm on SegTrack, we suggest you use our
cross-validation trained weights. Assuming that <DIR> is the directory where you
uncompressed the release archive, in the Matlab prompt type:

		cd <DIR>
    compile   #compiles all the nescessary mex files
    segmentSegTrack
    
Note: The release archive already contains precompiled mex files for 64-bit
Linux, so you only need to run the compile.m script if you are on another
platform (e.g. 32-bit Linux or Windows).

The outputs (optical flow, superpixels, segmentations, visuals) will be stored
in ‘<DIR>/Data/outputs/SegTrack/’. Note, that since the paper was submitted, we
have improved the algorithm by adding auto-calibration for some of the
parameters, so the segmentations produced will differ from the numbers reported
in [1]. 


Changelog
--------------------------------------------------------------------------------
v1.2: The demo.m file now resizes the input images by default. The file 
      demo_resizeFrame.m that showcased that feature before has been removed.

v1.1: Added script (segmentSegTrack.m) to run algorithm on SegTrack.
			Added script (demo_resizeFrames.m) to show how to run the code with online
			frame resizing in order to speed up (mainly) the optical flow computation.
			
v1.0.1: Removed option to compile the wrapper for [3]. The authors of [3] only
      provide binaries for 64-bit Linux, and a precompiled wrapper for that
      platform is already provided in the release.

v1.0: First public release



References
--------------------------------------------------------------------------------
[1] Anestis Papazoglou and Vittorio Ferrari
    Fast object segmentation in unconstraint video,
    ICCV 2013, Sydney, Australia

[2] Thomas Brox and Jitendra Malik
    Large Displacement Optical Flow: Descriptor Matching in Variational Motion Estimation,
    PAMI March 2011

[3] Narayanan Sundaram, Thomas Brox and Kurt Keutzer
    Dense Point Trajectories by GPU-accelerated Large Displacement Optical Flow,
    ECCV 2010, Crete, Greece

[4] Alex Levinshtein, Adrian Stere, Kiriakos N. Kutulakos, David J. Fleet, Sven J. Dickinson and Kaleem Siddiqi
    TurboPixels: Fast Superpixels Using Geometric Flows,
    PAMI 2009

[5] Radhakrishna Achanta, Appu Shaji, Kevin Smith, Aurelien Lucchi, Pascal Fua, and Sabine Süsstrunk
    SLIC Superpixels Compared to State-of-the-art Superpixel Methods,
    PAMI May 2012

[6] Pushmeet Kohli and Philip H.S. Torr
    Efficiently Solving Dynamic Markov Random Fields Using Graph Cuts,
    ICCV 2005, Beijing, China

[7] Alessandro Prest, Christian Leistner, Javier Civera, Cordelia Schmid and Vittorio Ferrari
    Learning Object Class Detectors from Weakly Annotated Video,
    CVPR 2012, Providence


