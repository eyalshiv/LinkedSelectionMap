

rmpath(genpath('E:/code/matlab'));
addpath(genpath('E:/GitHub/LinkedSelectionMap/matlab'));
addpath('E:/GitHub/LinkedSelectionMap/example/matlab');

global MLParamsStruct;

gl_initMLParamsStruct();


%% load basic information: chromosomes data, genetic maps


runLS_filenames_example;


% for c=1:length(genmap_files)
%   genmap{c}.file = genmap_files{c};
%   genmap{c}.name = genmap_token;
%   [genmap{c}.pos, genmap{c}.c, genmap{c}.R]  = textread( genmap{c}.file, '%d\t%f\t%f', 'commentstyle', 'shell' );
% end



cfg_inf  = LS_DefaultConfiguration( infcfg_file );

cfg_inf.inf.SWanno2param_mapping  = [1 2 3];
cfg_inf.inf.BSanno2param_mapping  = [1 2 3 4];
cfg_inf.inf.chromosomes           = cfg_inf.inf.chromosomes(1);
cfg_inf.inf.predef_idx_train      = cfg_inf.inf.predef_idx_train(1);
cfg_inf.chr_id                    = cfg_inf.chr_id(1);
cfg_inf.predef_idx_test           = cfg_inf.predef_idx_test(1);
cfg_inf.chromosomes               = cfg_inf.chromosomes(1);
cfg_inf.output.chromosomes        = cfg_inf.output.chromosomes(1);
cfg_inf.GEs.CalcSW.FE_grid        = 10.^-[1.5:1:4.5];
cfg_inf.GEs.CalcBS.FE_grid        = 10.^-[1.5:1:4.5];


%% build Grid Elements (GEs) for BS and SW 

struct2file( files_buildGE, files_buildGE_file );

cfg_inf.GEs.CalcSW.skip_generate_maps = 0;
cfg_inf.GEs.CalcBS.skip_generate_maps = 0;
struct2file( cfg_inf, infcfg_file );

% this can be skipped because 'LS_InferModel_f' builds the GEs if configured to
LS_PrecalcGridElements( files_buildGE_file, infcfg_file );



%% run inference

% since we built the GEs, don't do it again - just load them
cfg_inf.GEs.CalcSW.skip_generate_maps = 1;
cfg_inf.GEs.CalcBS.skip_generate_maps = 1;
struct2file( cfg_inf, infcfg_file );

struct2file( files_invar, files_invar_file );

struct2file( files_masks, files_masks_file );

[calco, mapo, outputo] = LS_InferModel_f( 'LS_BS4SW3_', files_invar_file, files_buildGE_file, files_masks_file, infcfg_file );
 


  
  






