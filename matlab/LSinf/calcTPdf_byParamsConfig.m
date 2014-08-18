
function [v_t, w_t, w_t_rel, u_del_max, u_del, u_del_rel, nsites_del] = calcTPdf_byParamsConfig( t_params, BSbase_params, config )
% function [v_t, w_t, u_del, nsites_del] = calcTPdf_byParamsConfig( t_params, BSbase_params, config )

global MLParamsStruct;



v_t = [];
w_t = [];
w_t_rel = [];
u_del = [];

% minimal_log10_t = -10;
minimal_log10_t = -10 + 0.001;

u_del_max = config.u_del_max;

nsites_del = [];
for a=1:length(BSbase_params)
  nsites_del(a) = BSbase_params{a}.nsites;
end

% switch config.t_pdf

% a mixture of point masses, for a single annotation:
% the first parameter marks the annotation used and the rest the mixture weights
%   case 'atoms'
for a=1:length(BSbase_params)
  
  aa = config.BSanno2param_mapping(a);
  
  v_t{a} = BSbase_params{a}.t;
  nvectors = length(v_t{a});
  assert( nvectors <= MLParamsStruct.bsparam_masses ,'precalculated base contains more vectors than possible degrees of freedom' );
  
  cur_weights = t_params((aa-1)*MLParamsStruct.bsparam_annolen+[1:nvectors]);
  irangeL = find(cur_weights <= minimal_log10_t);
  absw_t = 10.^cur_weights;
  absw_t(irangeL) = 0;
  if sum(absw_t)>0
    w_t{a} = absw_t/sum(absw_t);
  else
    w_t{a} = zeros(size(absw_t));
  end
  u_del(a)     = sum(absw_t);
  u_del_rel(a) = u_del(a) / u_del_max;
  
  w_t_rel{a} = absw_t ./ BSbase_params{a}.u;
end

% end

end

