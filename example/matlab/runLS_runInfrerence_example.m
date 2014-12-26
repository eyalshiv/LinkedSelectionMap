

% Making sure only the relevant GIT path is used
rmpath(genpath('E:/code/matlab'));
addpath(genpath('E:/GitHub/LinkedSelectionMap/matlab'));
% Adding to the path this example's local directory
addpath('E:/GitHub/LinkedSelectionMap/example/matlab');

% This global struct holds the description of the structure of the
% inference parameters vector (i.e. the mapping of each parameter to
% indices at the 106-parameters vector)
global MLParamsStruct;

gl_initMLParamsStruct();


%% Load basic information: chromosomes data, genetic maps

% Load the struct containing all input/output file names
runLS_filenames_example;

% Create a default configuration struct (note that this struct has all parameters required to run the inference besides the input/output file names) 
cfg_inf  = LS_DefaultConfiguration( infcfg_file );

% Use 3 annotations for CS
cfg_inf.inf.SWanno2param_mapping  = [1 2 3];
% Use 4 annotations for BS
cfg_inf.inf.BSanno2param_mapping  = [1 2 3 4];
% Run the inference on one chromosome only (2L)
cfg_inf.inf.chromosomes           = cfg_inf.inf.chromosomes(1);
cfg_inf.inf.predef_idx_train      = cfg_inf.inf.predef_idx_train(1);
cfg_inf.chr_id                    = cfg_inf.chr_id(1);
cfg_inf.predef_idx_test           = cfg_inf.predef_idx_test(1);
cfg_inf.chromosomes               = cfg_inf.chromosomes(1);
cfg_inf.output.chromosomes        = cfg_inf.output.chromosomes(1);
% Use a 4 masses distribution of fitness  effects for each annotation, both
% for CS and BS
cfg_inf.GEs.CalcSW.FE_grid        = 10.^-[1.5:1:4.5];
cfg_inf.GEs.CalcBS.FE_grid        = 10.^-[1.5:1:4.5];


%% Build Grid Elements (GEs) for BS and SW 

% Save input/output configuration struct to a file
struct2file( files_buildGE, files_buildGE_file );

% Force creation of grid elements
cfg_inf.GEs.CalcSW.skip_generate_maps = 0;
cfg_inf.GEs.CalcBS.skip_generate_maps = 0;
% Save inference configuration struct to a file
struct2file( cfg_inf, infcfg_file );

% Pre-calculating the grid elements
% (This can be skipped because 'LS_InferModel_f' automatically builds the GEs if
% configured to. It appears here to exemplify how to separate the creation of the GEs from the inference. )
LS_PrecalcGridElements( files_buildGE_file, infcfg_file );



%% Run inference

% Since we built the GEs, don't do it again - reconfigure to just load them
cfg_inf.GEs.CalcSW.skip_generate_maps = 1;
cfg_inf.GEs.CalcBS.skip_generate_maps = 1;
struct2file( cfg_inf, infcfg_file );

% Save variation input file names struct to file
struct2file( files_invar, files_invar_file );
% Save mask file names struct to file
struct2file( files_masks, files_masks_file );

% Finally, run the inference given all configuration files
[calco, mapo, outputo] = LS_InferModel_f( 'LS_BS4SW3_', files_invar_file, files_buildGE_file, files_masks_file, infcfg_file );
 


  
  






