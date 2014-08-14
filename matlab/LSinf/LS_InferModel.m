
function calc  = LS_InferModel( outfile_pref, fdata, GEs, cfg_inf, masks )

global MLParamsStruct;


opts_Alg{1} = optimset('Algorithm','active-set',     'Display','iter');
opts_Alg{2} = optimset('Algorithm','interior-point', 'Display','iter');
opts_Alg{3} = optimset('Algorithm','sqp',            'Display','iter');


%% inits

if isempty( cfg_inf )
  te      = LS_DefaultConfiguration( [outfile_pref 'infcfg.txt'], cfg_inf );
  cfg_inf = te.inf;
end



%% prepare data subset




%% 

% genome-wide average basic data statistics, e.g. mean heterozygosity, probability of a segregating site, synonymous divergence.
% NOTE: the statistics are calculated across the whole genome and not just
% for the specific codons used for the inference. I think for our purposes
% here it's good enough, but this may have to be changed.
gwEstats = genomewideStatisticsLight( fdata(cfg_inf.chromosomes) );


% Initializing parameters, either based on data assuming no selection, or from external parameters

tau0 = (gwEstats.Het/gwEstats.MutProx)^-1;

bnd_t_dist.fixed = cfg_inf.fixed_params;

bnd_t_dist.rangeL(1:MLParamsStruct.length) = 0;
bnd_t_dist.rangeL(MLParamsStruct.tau_pos)  = tau0 / 10;
bnd_t_dist.rangeL(MLParamsStruct.swparam_offset+[1:MLParamsStruct.swparam_annotations*(MLParamsStruct.swparam_masses+1)]) = 0;
bnd_t_dist.rangeL(MLParamsStruct.bsparam_offset+[1:MLParamsStruct.bsparam_annotations*(MLParamsStruct.bsparam_masses+1)]) = MLParamsStruct.minimal_log10_t;
bnd_t_dist.rangeL(MLParamsStruct.bsparam_imaxu) = log10(cfg_inf.u_del_max);
bnd_t_dist.rangeL(MLParamsStruct.swparam_imaxu) = 1;

bnd_t_dist.rangeH(MLParamsStruct.tau_pos)  = tau0 * 10;
bnd_t_dist.rangeH(MLParamsStruct.swparam_offset+[1:MLParamsStruct.swparam_annotations*(MLParamsStruct.swparam_masses+1)]) = 1;
bnd_t_dist.rangeH(MLParamsStruct.bsparam_offset+[1:MLParamsStruct.bsparam_annotations*(MLParamsStruct.bsparam_masses+1)]) = log10(10^-6); % I know its -6, just to make things clear
bnd_t_dist.rangeH(MLParamsStruct.bsparam_imaxu) = log10(cfg_inf.u_del_max);
bnd_t_dist.rangeH(MLParamsStruct.swparam_imaxu) = 1;

if isempty(cfg_inf.init_params)
  bnd_t_dist.init                         = bnd_t_dist.rangeL;
  bnd_t_dist.init(MLParamsStruct.tau_pos) = tau0;
else
  bnd_t_dist.init = min( max( cfg_inf.init_params, bnd_t_dist.rangeL ), bnd_t_dist.rangeH );
end


% for non-fixed selection weight params, make sure they are not null, since initial null values sometimes stuck the optimization
for k=1:MLParamsStruct.bsparam_annotations
  pidx = intersect( MLParamsStruct.bsparam_imasses(k)+[0:MLParamsStruct.bsparam_masses-1], find(~bnd_t_dist.fixed) );
  bnd_t_dist.init(pidx) = max( -9, bnd_t_dist.init(pidx) );
end
for k=1:MLParamsStruct.swparam_annotations
  pidx = intersect( MLParamsStruct.swparam_imasses(k)+[0:MLParamsStruct.swparam_masses-1], find(~bnd_t_dist.fixed) );
  bnd_t_dist.init(pidx) = max( 10^-8, bnd_t_dist.init(pidx) );
end


% fix&null all weights that do not correspond to a vector in GEs (in case we have less than 11 point masses per some annotation)
for k=1:MLParamsStruct.bsparam_annotations
  pidx = MLParamsStruct.bsparam_imasses(k)+[length(GEs.BSbase{1,k}.cfg.FE_grid):MLParamsStruct.bsparam_masses-1];
  bnd_t_dist.fixed(pidx) = 1;
  bnd_t_dist.init(pidx)  = MLParamsStruct.minimal_log10_t;
end
for k=1:MLParamsStruct.swparam_annotations
  pidx = MLParamsStruct.swparam_imasses(k)+[length(GEs.SWbase{1,k}.cfg.FE_grid):MLParamsStruct.swparam_masses-1];
  bnd_t_dist.fixed(pidx) = 1;
  bnd_t_dist.init(pidx)  = 0;
end

%% 
 
% % These are 'just-in-case' old lines, to remind me past mistakes...
%   cfg_inf.rec_th_L = anals_obs{r}.config.rec_th_L;
%   cfg_inf.rec_pol_th_L = cfg_inf.rec_th_L;


  calc         = inferNFSW( GEs, fdata, cfg_inf, bnd_t_dist, masks, []);

  struct2file( calc, [outfile_pref 'vparams.txt'] );

  bigo = 7;

end

