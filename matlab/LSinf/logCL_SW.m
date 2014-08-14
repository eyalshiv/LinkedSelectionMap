
function [neg_log_P, neg_log_samples, stats, vlogL] = logCL_SW(preCalc, full_params, config, variables, ivariables)

global MLParamsStruct;

minRed = 0.001;

calc_stats = (nargout > 2);


if nargin < 5
  variables  = [];
  ivariables = [];
end

params             = full_params;
params(ivariables) = variables;


iC = config.chromosomes;


TAU_div     = params(3);

if sum(isnan(params))
  neg_log_P = 100000000000000000000;
  return;
end

neg_log_P = 0;
neg_log_samples = [0 0 0];

for c=iC
  
  DivRedPred = composeLSMapFromElements( 1, preCalc{c}.gSWj, preCalc{c}.gBSj, params, preCalc{c}.SWparams, preCalc{c}.BSparams, preCalc{c}.EgMutDiv, config );
  
  theta0   = preCalc{c}.gMutDiv  / TAU_div; 
  
  pi0 = max( theta0./(1+theta0), 10^-5);
  
  pii        = min( pi0 .* max( DivRedPred.Red, minRed ), 0.99 );
  
  vlogL{c}      = double(preCalc{c}.samplesHom).*log(1-pii) + double(preCalc{c}.samplesHet).*log(pii);
  cneg_log_P{c}  = -sum(vlogL{c});
  
  cneg_log_samples{c} = [preCalc{c}.samplesHomSum+preCalc{c}.samplesHetSum  preCalc{c}.samplesHomSum  preCalc{c}.samplesHetSum];
  csites{c} = length(preCalc{c}.samplesHet);
  
  if nargout>2
    SumMutProx(c) = preCalc{c}.SumMutProx;
    nsamples(c) = preCalc{c}.samplesSum;
    nhetpairs(c)= preCalc{c}.samplesHetSum;
  end
end


for c=iC
  neg_log_samples = neg_log_samples + cneg_log_samples{c};
  neg_log_P       = neg_log_P       + cneg_log_P{c};
end

neg_log_P = neg_log_P/neg_log_samples(1) * 10^5;

% if isnan(neg_log_P) | isinf(neg_log_P)
%   bigo = 7;
% end
% [params(3) params(14:18)  neg_log_P/10000]
% [params(1) params(3) params(7:end) neg_log_P/1000000] %;params(1) params(3) w_s neg_log_P/1000000
% [params(1:3) neg_log_P/1000000]

if nargout>2
  stats.sites       = csites;
  stats.samplepairs = cneg_log_samples;
  stats.nlogLH      = cneg_log_P;
  stats.EMutProx    = sum(SumMutProx)  / sum(nsamples);
  stats.EHet        = sum(nhetpairs) / sum(nsamples);
end


bigo = 7;

end
