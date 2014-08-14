
function params = integrateParams_Inference_N_GEs( Inf_params, Inf_config, BSbase_params, SWbase_params, isX, iFE_grid_X )

global MLParamsStruct;


params.tau_div       = Inf_params(MLParamsStruct.tau_pos);

[ params.BS.v_t,...
  params.BS.w_t,...
  params.BS.rel_u_t,...
  params.BS.u_del,...
  params.BS.nsites_del] = calcTPdf_byParamsConfig(   Inf_params(MLParamsStruct.bsparam_offset+[1:MLParamsStruct.bsparam_annolen*MLParamsStruct.bsparam_annotations]), BSbase_params, Inf_config );

for a=1:length(params.BS.v_t)
  params.BS.U_del(a)    = params.BS.u_del(a)*params.BS.nsites_del(a);
  params.BS.Et(a)       =       sum(params.BS.w_t{a}.*params.BS.v_t{a} );
  params.BS.St(a)       = sqrt( sum(params.BS.w_t{a}.*params.BS.v_t{a}.^2) - params.BS.Et(a)^2 );
end



[ params.SW.v_s,...
  params.SW.coalrate_s] = calcSPdf_byParamsConfigNF( Inf_params(MLParamsStruct.swparam_offset+[1:MLParamsStruct.swparam_annolen*MLParamsStruct.swparam_annotations]), SWbase_params, Inf_config );
if ~isempty(SWbase_params)
  params.SW.Ne0               = SWbase_params{1}.Ne0;
  params.SW.remote2closeDiv   = SWbase_params{1}.remote2closeDiv;
  params.SW.tau_subs          = params.tau_div / params.SW.remote2closeDiv;
end
for a=1:length(params.SW.v_s)
  params.SW.alpha_over_tauS(a) = sum( params.SW.coalrate_s{a} );
  if params.SW.alpha_over_tauS(a) > 0
    params.SW.w_s{a} = params.SW.coalrate_s{a} / params.SW.alpha_over_tauS(a);
  else
    params.SW.w_s{a} = params.SW.coalrate_s{a};
  end
  
  params.SW.alpha(a)    = params.SW.alpha_over_tauS(a) * params.SW.tau_subs;
  
  params.SW.Es0(a)      =       sum(params.SW.alpha(a)*params.SW.w_s{a}.*params.SW.v_s{a});
  params.SW.Ss0(a)      = sqrt( sum(params.SW.alpha(a)*params.SW.w_s{a}.*params.SW.v_s{a}.^2) - params.SW.Es0(a)^2 );
  params.SW.Es(a)       =       sum(params.SW.w_s{a}.*params.SW.v_s{a} );
  params.SW.Ss(a)       = sqrt( sum(params.SW.w_s{a}.*params.SW.v_s{a}.^2) - params.SW.Es(a)^2 );  
end




if nargin < 4
  isX = 0;
end
if nargin < 5
  iFE_grid_X = [1:MLParamsStruct.bsparam_masses];
end

if isX
  
  for a=1:length(params.SW.alpha_over_tauS)
%     params.SW.w_s_X{a} = zeros(size(params.SW.w_s{a}));
    grid_s = params.SW.v_s{a}(iFE_grid_X(end:-1:1));
    ws     = params.SW.w_s{a}(iFE_grid_X(end:-1:1));
    [wsX, vsX_log10] = resampleDiscreteDistribution( log10( 4/3 * grid_s ), ws, log10( grid_s ),  log10( [grid_s(1)/10 grid_s(end)*10] ) );
    %   [wsX, vsX_log10] = resampleDiscreteDistribution(        4/3 * grid_s  , ws,        grid_s  ,         [grid_s(1)/10 grid_s(end)*10]   );
    params.SW.w_s_X{a} = wsX(end:-1:1);
    params.SW.v_s_X{a} = 10.^vsX_log10(end:-1:1);
  end
  
  for a=1:length(params.BS.u_del)
%     params.BS.w_t_X{a} = zeros(size(params.BS.w_t{a}));
    grid_t = params.BS.v_t{a}(iFE_grid_X(end:-1:1));
    wt     = params.BS.w_t{a}(iFE_grid_X(end:-1:1));
    [wtX, vtX_log10] = resampleDiscreteDistribution( log10( 4/3 * grid_t ), wt, log10( grid_t ),  log10( [grid_t(1)/10 grid_t(end)*10] ) );
    params.BS.w_t_X{a} = wtX(end:-1:1);
    params.BS.v_t_X{a} = 10.^vtX_log10(end:-1:1);
  end
  
end





end

