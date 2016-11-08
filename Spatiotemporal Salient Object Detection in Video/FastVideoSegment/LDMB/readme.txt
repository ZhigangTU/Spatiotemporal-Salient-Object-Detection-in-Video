###################################################################
#                                                                 #
#    Motion Boundaries Detection v1.0                             #
#    Philippe Weinzaepfel (philippe.weinzaepfel@inria.fr)         #
#                                                                 #
###################################################################

1. Introduction.

Code for detecting motion boundaries built upon structured random forest proposed in [Structured Forest for Fast Edge Detection, Dollar and Zitnick, ICCV'13].

If you use our code, please cite:

@inproceedings{weinzaepfel:hal-01142653,
  TITLE = {{Learning to Detect Motion Boundaries}},
  AUTHOR = {Weinzaepfel, Philippe and Revaud, Jerome and Harchaoui, Zaid and Schmid, Cordelia},
  BOOKTITLE = {{CVPR 2015 - IEEE Conference on Computer Vision \& Pattern Recognition}},
  YEAR = {2015},
}

###################################################################

2. License.

This code is published under the MSR-LA Full Rights License.
Please read license.txt for more info.

###################################################################

3. Installation.

a) This code is written for the Matlab interpreter (tested with versions R2014a-2014b) and requires the Matlab Image Processing Toolbox. 

b) Additionally, Piotr's Matlab Toolbox (version 3.26 or later) written by Piotr Dollar is also required. It can be downloaded at:
 http://vision.ucsd.edu/~pdollar/toolbox/doc/index.html.

c) Next, please compile mex code from within Matlab:
  mex private/edgesDetectMex.cpp -outdir private [OMPPARAMS]
  mex private/image_warping_error.cpp -outdir private [OMPPARAMS]
Here [OMPPARAMS] are parameters for OpenMP and are OS and compiler dependent.
  Windows:  [OMPPARAMS] = '-DUSEOMP' 'OPTIMFLAGS="$OPTIMFLAGS' '/openmp"'
  Linux: [OMPPARAMS] = '-DUSEOMP' CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
To compile without OpenMP simply omit [OMPPARAMS]; note that code will be single threaded in this case.

d) Models based on classic+nl-fast optical flow estimation are provided. The code for computing the flow can be downloaded at:
 http://cs.brown.edu/~dqsun/research/software.html

###################################################################

4. Getting Started.

 - Make sure to carefully follow the installation instructions above.
 - Please see "demo.m"

###################################################################

5. History.

Version 1.0 (06/2015)

###################################################################
