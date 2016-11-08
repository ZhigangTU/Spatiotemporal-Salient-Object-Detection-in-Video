% The make utility for all the C and MEX code

function make(command)

if (nargin > 0 && strcmp(command,'clean'))
    delete('*.mexglx');
    delete('*.mexw32');
    delete('lsmlib/*.mexglx');
    delete('lsmlib/*.mexw32');
    return;
end
mex DT.c
mex height_function_der.c
mex height_function_grad.c
mex local_min.c
mex zero_crossing.c
mex -lm get_full_speed.c
mex corrDn.c wrap.c convolve.c edges.c
mex upConv.c wrap.c convolve.c edges.c

cd lsmlib
mex computeDistanceFunction2d.c FMM_Core.cpp FMM_Heap.cpp lsm_FMM_field_extension2d.cpp
mex computeExtensionFields2d.c FMM_Core.cpp FMM_Heap.cpp lsm_FMM_field_extension2d.cpp
mex doHomotopicThinning.c FMM_Core.cpp FMM_Heap.cpp lsm_FMM_field_extension2d.cpp
cd ..