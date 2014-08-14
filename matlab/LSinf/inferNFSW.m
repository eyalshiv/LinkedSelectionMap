
function [calc, preLHstats]  = inferNFSW( GEs, fdata, cfg_inf, bounds, masks, preLHstats)


% cfg_inf.fmincon_retries = 0;


calc.init_params = bounds.init;
calc.rangeL = bounds.rangeL;
calc.rangeH = bounds.rangeH;
calc.rangeL(bounds.fixed==1) = bounds.init(bounds.fixed==1);
calc.rangeH(bounds.fixed==1) = bounds.init(bounds.fixed==1);

calc.DFs = sum( calc.rangeL<calc.rangeH );


iC = cfg_inf.chromosomes;



if isfield(cfg_inf,'constraint_u') && (cfg_inf.constraint_u == 1)
  f_constraint = @constraint_u_sum;
elseif isfield(cfg_inf,'constraint_u') && (cfg_inf.constraint_u == 2)
  f_constraint = @constraint_u_sum_ineq;
else
  f_constraint = [];
end

if ~isfield(cfg_inf, 'opt_type') | isempty(cfg_inf.opt_type)
  cfg_inf.opt_type = {'active-set', 'interior-point', 'sqp'};
end

for i=1:length(cfg_inf.opt_type)
  opts{i} = optimset('Algorithm',cfg_inf.opt_type{i},     'Display','iter');
end

if isempty(preLHstats)
  for c=iC
%     if isfield(cfg_inf, 'predef_idx_train') & ~isempty(cfg_inf.predef_idx_train)
    if exist('masks') & ~isempty(masks)
      predef_idx = masks{c};
    else
      predef_idx = [];
    end
    
    if ~isfield(GEs, 'BSbase')
      cur_BSbase = {[]};
    else
      cur_BSbase = GEs.BSbase(c,:);
    end
    if ~isfield(GEs, 'SWbase')
      cur_SWbase = {[]};
    else
      cur_SWbase = GEs.SWbase(c,:);
      cur_gFocGrid = GEs.gFocGrid{c};
    end
    
    preLHstats{c}    = calcPreLHStatsSw( cur_SWbase, cur_gFocGrid, cur_BSbase, fdata{c}, cfg_inf, predef_idx );
    %     calc.EgMutDiv(c) = preLHstats{c}.EgMutDiv;
  end
end

calc.init_LH = logCL_SW(preLHstats, calc.init_params, cfg_inf);
calc.best_iter = -1;
LH_best = calc.init_LH;
% params_temp  = calc.init_params;
% LH_temp      = calc.init_LH;

opts_default = optimoptions('fmincon');

for k=1:length(opts)
  
  if calc.best_iter==-1
    calc.params(k,:)  = calc.init_params;
  else
    calc.params(k,:)  = calc.params(calc.best_iter,:);
  end
  
  retries = cfg_inf.fmincon_retries;
  calc.exitflag(k)    = 0;
  opts{k}.MaxIter     = opts_default.MaxIter;
  opts{k}.MaxFunEvals = 100*size(calc.params,2); %opts_default.MaxFunEvals;
  
  while ( retries >= 0 ) & ( calc.exitflag(k) == 0 )
    
    full_params = calc.params(k,:);
    ivariables  = find(bounds.fixed==0);
  
    if strcmp(cfg_inf.opt_type{k}, 'interior-point')
      full_params(ivariables) = epsilonFromBoundary( full_params(ivariables), calc.rangeL(ivariables), calc.rangeH(ivariables), 10^-9 );
    end
 
    
    if isfield(cfg_inf,'constraint_u') && (cfg_inf.constraint_u == 1)
      [tvariables, calc.nLogLH(k), calc.exitflag(k), calc.opt{k}] = fmincon(@(variables)logCL_SW(preLHstats, full_params, cfg_inf, variables, ivariables), full_params(ivariables), [],[],[],[], calc.rangeL(ivariables), calc.rangeH(ivariables), @(variables)constraint_u_sum(full_params,variables,ivariables), opts{k});
    elseif isfield(cfg_inf,'constraint_u') && (cfg_inf.constraint_u == 2)
      [tvariables, calc.nLogLH(k), calc.exitflag(k), calc.opt{k}] = fmincon(@(variables)logCL_SW(preLHstats, full_params, cfg_inf, variables, ivariables), full_params(ivariables), [],[],[],[], calc.rangeL(ivariables), calc.rangeH(ivariables), @(variables)constraint_u_sum_ineq(full_params,variables,ivariables), opts{k});
    else
      [tvariables, calc.nLogLH(k), calc.exitflag(k), calc.opt{k}] = fmincon(@(variables)logCL_SW(preLHstats, full_params, cfg_inf, variables, ivariables), full_params(ivariables), [],[],[],[], calc.rangeL(ivariables), calc.rangeH(ivariables), [], opts{k});
    end
    calc.params(k,ivariables) = tvariables;
    
    if calc.exitflag(k) == 0
      opts{k}.MaxFunEvals = 2*opts{k}.MaxFunEvals; %[num2str(2*str2num(opts{k}.MaxFunEvals(1:find(opts{k}.MaxFunEvals=='*')-1))) '*numberOfVariables'];
      opts{k}.MaxIter     = 2*opts{k}.MaxIter;
    end
    retries = retries-1;
  end
  
  if ( LH_best >= calc.nLogLH(k) ) & ( calc.exitflag(k) >= 0 )
    calc.best_iter = k;
    LH_best = calc.nLogLH(k);
  end
  %   if LH_temp >= calc.nLogLH(k)  &  calc.exitflag(k) > 0
  % %   if LH_temp >= calc.nLogLH(k)
  %     params_temp = calc.params{k};
  %     LH_temp     = calc.nLogLH(k);
  %     calc.best_iter = k;
  %   end
  
end

calc.nLogLH = calc.nLogLH / 10^5;

if calc.best_iter~=-1
  calc.best_iter = find( calc.nLogLH == min(calc.nLogLH), 1, 'last' );
  [~, calc.samples, calc.LHFuncStats] = logCL_SW(preLHstats, calc.params(calc.best_iter,:), cfg_inf);
  calc.logLH = -calc.nLogLH( calc.best_iter );
else
  calc.samples = [0 0 0];
  calc.logLH = 0;
end



% if isfield(cfg_inf, 'predef_idx_train')
%   % this is a patch in order to save space - theoreticaly we need to track
%   % the exact list of codons used for the inference, but we settle with
%   % keeping this field only at the 'anals_obs' struct for all corresponding
%   % inference runs
%   calc.config = rmfield(cfg_inf, 'predef_idx_train');
% else
%   bigo = 7;
% end
% % if isfield(calc.config, 'predef_idx_test')
% %   calc.config = rmfield(calc.config, 'predef_idx_test');
% % end
% % if isfield(calc.config, 'pos')
% %   calc.config = rmfield(calc.config, 'pos');
% % end



calc.config.inf = cfg_inf;
calc.config.GEs = GEs.cfg;


end

