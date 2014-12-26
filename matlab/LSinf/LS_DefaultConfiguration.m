
function ncfg   = LS_DefaultConfiguration( outfile, cfg )

% This function creates a new full inference configuration struct or
% updates/completes a given one.

% inputs:
% outfile    - output configuration file name
% cfg        - existing configuration struct
% outputs:
% ncfg       - new/updated configuration struct


global MLParamsStruct;
gl_initMLParamsStruct();


% Initialize all configuration fields (for both data processing, model inferece and model evaluation).
% If 'cfg' has some of the fields, only the missing ones are initialized.

if nargin < 2
  ncfg = [];
else
  ncfg = cfg;
end

if isstr(ncfg)
  ncfg = file2struct(ncfg);
end


if ~isfield(ncfg,     'inf'),                   ncfg.inf= []; ,end
if ~isfield(ncfg.inf, 'optim_algs'),            ncfg.inf.optim_algs = [1 2 3]; ,end
if ~isfield(ncfg.inf, 'fmincon_retries '),      ncfg.inf.fmincon_retries = 0; ,end
if ~isfield(ncfg.inf, 'chromosomes'),           ncfg.inf.chromosomes = [1 2 3 4]; ,end
if ~isfield(ncfg.inf, 'constraint_u'),          ncfg.inf.constraint_u = 2; ,end
if ~isfield(ncfg.inf, 'use_fake_subs'),         ncfg.inf.use_fake_subs     = 1;              ,end
if ~isfield(ncfg.inf, 'SWanno2param_mapping'),  ncfg.inf.SWanno2param_mapping = [1 2 3 4]; ,end
if ~isfield(ncfg.inf, 'BSanno2param_mapping'),  ncfg.inf.BSanno2param_mapping = [1 2 3 4]; ,end
if ~isfield(ncfg.inf, 'tau_sweeps'),            ncfg.inf.tau_sweeps = -1; ,end % when this is -1, calculate by scaling tau_mutprox, using the parameters 'ratio_tau_mutprox_tau_sweeps'
if ~isfield(ncfg.inf, 'ratio_tau_mutprox_tau_sweeps'), ncfg.inf.ratio_tau_mutprox_tau_sweeps = 3.25; ,end % default = 1 ?
if ~isfield(ncfg.inf, 'u_del_max'),             ncfg.inf.u_del_max = 6.8*10^-9; ,end
if ~isfield(ncfg.inf, 'remote2closeDiv'),       ncfg.inf.remote2closeDiv = 3.2583; ,end
if ~isfield(ncfg.inf, 'predef_idx_train'),            ncfg.inf.predef_idx_train = {[],[],[],[]}; ,end
if ~isfield(ncfg.inf, 'init_params'),        ncfg.inf.init_params = []; ,end
% if ~isfield(ncfg, 'iannoSW'),        ncfg.iannoSW = [1 2 3 4]; ,end
% if ~isfield(ncfg, 'iannoBS'),        ncfg.iannoBS = [1 2 3 4]; ,end
if ~isfield(ncfg.inf, 'fixed_params')
  ncfg.inf.fixed_params(1:MLParamsStruct.length) = 1;
  ncfg.inf.fixed_params(MLParamsStruct.tau_pos ) = 0;
  ncfg.inf.fixed_params([MLParamsStruct.bsparam_imasses(1) + [2 4 6 8 10]-1  MLParamsStruct.bsparam_imasses(2) + [2 4 6 8 10]-1  MLParamsStruct.bsparam_imasses(3) + [2 4 6 8 10]-1  MLParamsStruct.bsparam_imasses(4) + [2 4 6 8 10]-1]) =  0;
  ncfg.inf.fixed_params([MLParamsStruct.swparam_imasses(1) + [2 4 6 8 10]-1  MLParamsStruct.swparam_imasses(2) + [2 4 6 8 10]-1  MLParamsStruct.swparam_imasses(3) + [2 4 6 8 10]-1  MLParamsStruct.swparam_imasses(4) + [2 4 6 8 10]-1]) =  0;
end
if ~isfield(ncfg.inf, 'opt_type'), ncfg.inf.opt_type = {'active-set', 'interior-point', 'sqp'};   ,end


if ~isfield(ncfg.inf, 'chr_id'),               ncfg.chr_id = {'2L', '2R', '3L', '3R', 'X'}; ,end


% if ~isfield(ncfg.inf, 's_pdf'),        ncfg.s_pdf = 'atoms'; ,end
% if ~isfield(ncfg.inf, 't_pdf'),        ncfg.t_pdf = 'atoms'; ,end

if ~isfield(ncfg, 'predef_idx_test'),            ncfg.predef_idx_test = {[],[],[],[]}; ,end

if ~isfield(ncfg, 'chromosomes'),        ncfg.chromosomes = [1 2 3 4]; ,end

if ~isfield(ncfg, 'foc_rec_th_L'),        ncfg.foc_rec_th_L = 0.75; ,end
if ~isfield(ncfg, 'foc_rec_th_H'),        ncfg.foc_rec_th_H = Inf; ,end
if ~isfield(ncfg, 'rec_th_L'),        ncfg.rec_th_L = 0.75; ,end
if ~isfield(ncfg, 'rec_th_H'),        ncfg.rec_th_H = Inf; ,end
if ~isfield(ncfg, 'rec_pol_th_L'),        ncfg.rec_pol_th_L = 0.75; ,end
if ~isfield(ncfg, 'rec_pol_th_H'),        ncfg.rec_pol_th_H = Inf; ,end

if ~isfield(ncfg, 'gL_COL'),        ncfg.gL_COL = 0.0011; ,end
if ~isfield(ncfg, 'gdelta_COL'),        ncfg.gdelta_COL = 10^-6; ,end
if ~isfield(ncfg, 'sort_strand_dir'),        ncfg.sort_strand_dir = 1; ,end
if ~isfield(ncfg, 'halfway'),        ncfg.halfway = 0; ,end
if ~isfield(ncfg, 'collated_weights'),        ncfg.collated_weights = 0; ,end

%     if ~isfield(ncfg, 'div_smooth_idx'),        ncfg.div_smooth_idx = 5; ,end
if ~isfield(ncfg, 'scale_dS'),        ncfg.scale_dS = 0.4; ,end
if ~isfield(ncfg, 'use_paml_dS'),        ncfg.use_paml_dS = 1; ,end
if ~isfield(ncfg, 'min_paml_codons'),        ncfg.min_paml_codons = 50; ,end
if ~isfield(ncfg, 'smoothgenmap'),        ncfg.smoothgenmap = 1; ,end
if ~isfield(ncfg, 'smoothgenmap_win'),        ncfg.smoothgenmap_win = 10^-6; ,end

if ~isfield(ncfg, 'rec_spat_window'),        ncfg.rec_spat_window = 10^5; ,end



if ~isfield(ncfg, 'genmap_name'),        ncfg.genmap_name = 'Comeron'; ,end



if ~isfield(ncfg, 'output'),             ncfg.output = []; ,end

if ~isfield(ncfg.output, 'LSmap_res'),          ncfg.output.LSmap_res = 100; ,end
if ~isfield(ncfg.output, 'chromosomes'),          ncfg.output.chromosomes = [1 2 3 4]; ,end
if ~isfield(ncfg.output, 'spatial_resolution'),          ncfg.output.spatial_resolution = 1000; ,end



if ~isfield(ncfg,        'GEs'),          ncfg.GEs = []; ,end

if ~isfield(ncfg.GEs,        'CalcSW'),          ncfg.GEs.CalcSW = []; ,end
if ~isfield(ncfg.GEs.CalcSW, 'StopSum'),         ncfg.GEs.CalcSW.StopSum       = 1;              ,end
if ~isfield(ncfg.GEs.CalcSW, 'InterpMethod'),    ncfg.GEs.CalcSW.InterpMethod  = 'linear';       ,end
if ~isfield(ncfg.GEs.CalcSW, 'trap_aprx'),       ncfg.GEs.CalcSW.trap_aprx     = 'diffusion';    ,end  %DURRETT1;
if ~isfield(ncfg.GEs.CalcSW, 'gMaxDistScaled'),  ncfg.GEs.CalcSW.gMaxDistScaled= 1;              ,end  %scale the radius of maximal summation in proportion to S
if ~isfield(ncfg.GEs.CalcSW, 'gMaxDist'),        ncfg.GEs.CalcSW.gMaxDist      = 1;              ,end  %if no scaling, absolute distance in morgans, otherwise the proportionality constant between S and maximal radius (0.1 is reasonable given the rule-of-thumb r=0.1s, 1 is very safe)
if ~isfield(ncfg.GEs.CalcSW, 'FE_grid'),         ncfg.GEs.CalcSW.FE_grid = 10.^-[1:0.5:6];              ,end
if ~isfield(ncfg.GEs.CalcSW, 'Ne0'),             ncfg.GEs.CalcSW.Ne0 = 1.1408*10^6;              ,end
if ~isfield(ncfg.GEs.CalcSW, 'skip_generate_maps'),  ncfg.GEs.CalcSW.skip_generate_maps = 1; ,end
if ~isfield(ncfg.GEs.CalcSW, 'min_mapdist_SW_spat_grid'), ncfg.GEs.CalcSW.min_mapdist_SW_spat_grid = 3*10^-8; ,end
 
if ~isfield(ncfg.GEs,        'CalcBS'),          ncfg.GEs.CalcBS = []; ,end
if ~isfield(ncfg.GEs.CalcBS, 'FE_grid'),             ncfg.GEs.CalcBS.FE_grid = 10.^-[1:0.5:6];              ,end
if ~isfield(ncfg.GEs.CalcBS, 'u_del'),               ncfg.GEs.CalcBS.u_del      = 3.5*10.^-9;              ,end
if ~isfield(ncfg.GEs.CalcBS, 'skip_generate_maps'),  ncfg.GEs.CalcBS.skip_generate_maps = 1; ,end
if ~isfield(ncfg.GEs.CalcBS, 'B_res'),               ncfg.GEs.CalcBS.B_res = 100; ,end
if ~isfield(ncfg.GEs.CalcBS, 'exe_file'),            ncfg.GEs.CalcBS.exe_file = 'E:\projects\sweeps\NeutFoc\BS_util\BSmap.exe'; ,end
% if ~isfield(ncfg.GEs.CalcBS, 't_dist_type'),       ncfg.GEs.CalcBS.t_dist_type = 'POINT'; ,end



if ~isempty(outfile)
  struct2file( ncfg, outfile );
end


end

