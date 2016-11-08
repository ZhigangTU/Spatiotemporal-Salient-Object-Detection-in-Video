
This code includes the detailed implementation of CVPR10 paper.

Ref:
Yigang Peng, Arvind Ganesh, John Wright, Wenli Xu, and Yi Ma. "RASL: Robust Alignment
by Sparse and Low-rank Decomposition for Linearly Correlated Images," the 23th IEEE Conference
on Computer Vision and Pattern Recognition (CVPR 10'), San Francisco, CA, USA, June 2010.


You may start an example by directly running 'go_*.m' for various data.

go_man.m         figure 1 in cvpr paper
go_dummy.m       figure 3 in cvpr paper
go_LFW.m         figure 4 in cvpr paper (one of 20 persons)
go_gore.m        figure 5 in cvpr paper
go_digits_3.m    figure 6 in cvpr paper
go_windows.m     figure 7 in cvpr paper

Have fun!

The code contains:
By default, the data are stored in the 'data' folder, and the results are stored in 'results' folder.
Each image file is also accompanied by a .mat file containing the initial coordinates of feature points. 
The 'RASL_toolbox' folder contains some necessary functions for file organization and image processing.

rasl_main.m         the main iteration loop for RASL
rasl_inner_apg.m    the APG algorithm for the inner convex optimization
rasl_inner_ialm.m   the inexact algorithm for the inner convex optimization
You may switch between the two algorithm by using the 'APGorALM_flag' in 'rasl_main.m'.
Note that, the inexact ALM algorithm is experiencally much faster than the APG algorithm. 
However, the APG algorithm is gurantted to converge, while the inexact algorithm not.
By default, we use the inexact algorithm for the sake of time.


For algorithm details, please read our technical reports, in which we explain
more details of RASL.

Any comments or questions, please contact 
Yigang Peng (pyg07@mails.tsinghua.edu.cn)
Arvind Ganesh (abalasu2@illinois.edu)


