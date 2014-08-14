
function [gSWj, gBSj, SWparams, BSparams, DivRedPred] = composePredictedDiversity2( SWbase, SWgpos, BSbase, params, config, pos, EgMutDiv ) %, remote2closeDiv
% function [DivRedPred, sumeps] = composePredictedDiversity2( preLH, params, config, gpos )

% calculating the predicted reduction in diversity due to linked selection
% (selective sweeps and background selection)


global MLParamsStruct;

pos = double(pos);

% prepare SW and BS grid elements for all spatial positions

gSWj = [];
if ~isempty(SWbase{1})
  for a=1:length(SWbase)
    for k=1:length(SWbase{a}.gSWj)
      idx = ~isnan(SWbase{a}.gSWj{k});
      gSWj{a}(k,:)   =                interp1(double(SWgpos.pos{k}(idx)),  SWbase{a}.gSWj{k}(idx),      pos(:)');
%     gSWj{a}(k,:)   =                interp1(SWgpos.gpos{k}(idx), SWbase{a}.gSWj{k}(idx),      gpos(:)');
      if config.use_fake_subs
        gSWj{a}(k,:) = gSWj{a}(k,:) + interp1(double(SWgpos.pos{k}(idx)),  SWbase{a}.gSWj_fake{k}(idx), pos(:)');
%       gSWj{a}(k,:) = gSWj{a}(k,:) + interp1(SWgpos.gpos{k}(idx), SWbase{a}.gSWj_fake{k}(idx), gpos(:)');
      end
      gSWj{a}(k,isnan(gSWj{a}(k,:))) = 0;
    end
  end
else
  gSWj{1} = ones(size(pos))';
end

gBSj = [];
if ~isempty(BSbase{1})
  for a=1:length(BSbase)
    for k=1:length(BSbase{a}.Bj)
      cX = cumsum([0 BSbase{a}.Lj{k}(1:end-1)]);
      cX = [0.75+cX; BSbase{a}.Lj{k}+0.25+cX];
      cB = double([BSbase{a}.Bj{k} BSbase{a}.Bj{k}])';
      gBSj{a}(k,:) = interp1(double(cX(:)), double(cB(:)), double(pos));
    end
  end
else
  gBSj{1} = ones(size(pos))';
end


% collect (selection) parameters of the grid elements

SWparams = [];
if ~isempty(SWbase{1})
  for a=1:length(SWbase)
    SWparams{a}.s    = SWbase{a}.cfg.FE_grid;
    SWparams{a}.Ne0  = SWbase{a}.cfg.Ne0; %SWbase{a}.params{k}.Ne0;
    SWparams{a}.remote2closeDiv = config.remote2closeDiv;
  end
else
  SWparams{1}.s = 1;
  SWparams{1}.Ne0 = 1;
  SWparams{1}.remote2closeDiv = 1;
end

BSparams = [];
if ~isempty(BSbase{1})
  for a=1:length(BSbase)
    BSparams{a}.u         = BSbase{a}.cfg.u_del;
    BSparams{a}.t         = BSbase{a}.cfg.FE_grid;
    BSparams{a}.nsites    = BSbase{a}.cfg.anno_len;
  end
else
  BSparams{1}.u         = 10^-10;
  BSparams{1}.t         = 1;
  BSparams{1}.nsites    = 1;
end



% construct BS, SW and combined linked selection maps from the grid elements and parameters

if ~isempty(params)

  DivRedPred = composeLSMapFromElements( 0, gSWj, gBSj, params, SWparams, BSparams, EgMutDiv, config );
  
  DivRedPred.pos = pos;

else
  
  DivRedPred = [];
  
end


end
