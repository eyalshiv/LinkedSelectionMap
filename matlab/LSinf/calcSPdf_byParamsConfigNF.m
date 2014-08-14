
function [v_s, coalrate_s] = calcSPdf_byParamsConfigNF(s_params, SWbase_params, config)

global MLParamsStruct;

v_s = [];
coalrate_s = [];

for a=1:length(SWbase_params)
  v_s{a} = SWbase_params{a}.s;
  coalrate_s{a} = zeros(size(v_s{a}));

  nvectors = length(v_s{a});
  assert( nvectors <= MLParamsStruct.swparam_masses ,'precalculated base contains more vectors than possible degrees of freedom' );
  
  % for each annotation, find the appropriate parameters set accodring to the mapping from the configutaion, iSWanno
  aa = config.SWanno2param_mapping(a);
  
%   switch config.s_pdf
%     case 'atoms'

      
      coalrate_s{a} = max(s_params((aa-1)*(MLParamsStruct.swparam_masses+1)+[1:nvectors]),0);
%   end
  
end

end
