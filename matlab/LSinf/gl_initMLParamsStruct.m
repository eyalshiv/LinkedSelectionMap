
function gl_initMLParamsStruct()

global MLParamsStruct;


MLParamsStruct.tau_pos = 3;
MLParamsStruct.swparam_offset = 10;
MLParamsStruct.swparam_masses = 11;
MLParamsStruct.swparam_annotations = 4;
MLParamsStruct.swparam_annolen = MLParamsStruct.swparam_masses + 1;
MLParamsStruct.bsparam_offset = MLParamsStruct.swparam_offset + MLParamsStruct.swparam_annolen*MLParamsStruct.swparam_annotations;
MLParamsStruct.bsparam_masses = 11;
MLParamsStruct.bsparam_annotations = 4;
MLParamsStruct.bsparam_annolen = MLParamsStruct.bsparam_masses + 1;
MLParamsStruct.length          = MLParamsStruct.bsparam_offset + (MLParamsStruct.bsparam_masses+1)*MLParamsStruct.bsparam_annotations;

MLParamsStruct.swparam_imasses = MLParamsStruct.swparam_offset  + [0:MLParamsStruct.swparam_annotations-1]*MLParamsStruct.swparam_annolen + 1;
MLParamsStruct.swparam_imaxu   = MLParamsStruct.swparam_imasses + MLParamsStruct.swparam_masses;
MLParamsStruct.bsparam_imasses = MLParamsStruct.bsparam_offset  + [0:MLParamsStruct.bsparam_annotations-1]*MLParamsStruct.bsparam_annolen + 1;
MLParamsStruct.bsparam_imaxu   = MLParamsStruct.bsparam_imasses + MLParamsStruct.bsparam_masses;


% a technical constant - SHOULD NOT BE HERE !!!
MLParamsStruct.minimal_log10_t = -10;


end