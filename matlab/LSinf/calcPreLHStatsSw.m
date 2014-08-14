
function preLHstats = calcPreLHStatsSw( SWbase, SWgpos, BSbase, data, config, weights, ccounts, pos_external, gpos_external ) %remote2closeDiv, 


% min_paml_codons = 100;
% scale_dS = 1; %0.4;


if ~exist('gpos_external')
  preLHstats.gpos_ext = [];
else
  preLHstats.gpos_ext = gpos_external;
end
if ~exist('pos_external')
  preLHstats.pos_ext = [];
else
  preLHstats.pos_ext = pos_external;
end
if ~exist('ccounts')
  ccounts = [];
end
if ~exist('weights')
  weights = ones(sizeof(data.Poly.pos));
end
if length(weights) ~= length(data.Poly.pos)
  error('predef_idx length different from Poly.pos length');
end


idx1p = data.Poly.idxS;
% idx1d = data.LongDiv.idxS;


% if config.rec_spat_window > 0
%   grat = spatAverageRecRat( data.genmap, data.Poly.pos, config.rec_spat_window );
% else
%   grat = data.Poly.grat;
% end
% mask_included_codons = (grat >= config.rec_th_L) & (grat <= config.rec_th_H);

% mask_included_codons  = ones(size(data.Poly.pos));
%   mask_included_codons(setdiff([1:length(mask_included_codons)], find(predef_idx>0))) = 0;

preLHstats.idxc = find(weights>0);

[te, preLHstats.idxc1p, iidx1p] = intersect(preLHstats.idxc, idx1p);
% [te, preLHstats.idxc1d, iidx1d] = intersect(preLHstats.idxc, idx1d);
preLHstats.idxc1p = preLHstats.idxc1p';
% preLHstats.idxc1d = preLHstats.idxc1d';

preLHstats.gpos = data.Poly.gpos(preLHstats.idxc);
% preLHstats.grat = data.Poly.grat(preLHstats.idxc);
preLHstats.pos  = data.Poly.pos( preLHstats.idxc);

samples                     = [weights(preLHstats.idxc).*double(data.Poly.sample(preLHstats.idxc)) zeros([length(preLHstats.idxc) 1]) ];
samples(preLHstats.idxc1p,:) = (weights(preLHstats.idxc(preLHstats.idxc1p))*[1 1]).*double(data.Poly.specS(iidx1p,:));

preLHstats.samplesHet = uint16(prod(samples, 2)); % m_by_n
preLHstats.samplesHom = uint16(sum(samples,2).*(sum(samples,2)-1)/2 - double(preLHstats.samplesHet)); % n_over_2 - m_by_n
preLHstats.samplesHetSum = sum(preLHstats.samplesHet);
preLHstats.samplesHomSum = sum(preLHstats.samplesHom);
preLHstats.samplesSum    = preLHstats.samplesHomSum+preLHstats.samplesHetSum;

% if isfield(config, 'use_paml_dS') & config.use_paml_dS
  dSidx = find( ~isnan(data.MutProx.dS) & data.MutProx.dS >= 0 ); %find(data.MutProx.paml_codons >= min_paml_codons);
  preLHstats.gMutDiv   = interp1(double(data.MutProx.pos(dSidx)), data.MutProx.dS(dSidx), double(preLHstats.pos)); %data.MutProx.dS(dSidx)*scale_dS
  preLHstats.gMutDiv(preLHstats.pos < data.MutProx.pos(dSidx(1  ))) = data.MutProx.dS(dSidx(1  )); %*scale_dS;
  preLHstats.gMutDiv(preLHstats.pos > data.MutProx.pos(dSidx(end))) = data.MutProx.dS(dSidx(end)); %*scale_dS;
% else
%   preLHstats.gMutDiv   = interp1(data.MutProx.pos, data.MutProx.div, double(preLHstats.pos));
%   preLHstats.gMutDiv(preLHstats.pos < min(data.MutProx.pos)) = data.MutProx.div(1);
%   preLHstats.gMutDiv(preLHstats.pos > max(data.MutProx.pos)) = data.MutProx.div(end);
% end
preLHstats.SumMutProx = sum(preLHstats.gMutDiv.*(double(preLHstats.samplesHet)+double(preLHstats.samplesHom)));
preLHstats.EgMutDiv = nanmean(preLHstats.gMutDiv);
preLHstats.EgHet    = preLHstats.samplesHetSum/(preLHstats.samplesHetSum+preLHstats.samplesHomSum);



if isempty(preLHstats.pos_ext)|isempty(preLHstats.gpos_ext)
  cur_pos  = double(data.Poly.pos(preLHstats.idxc));
%   cur_gpos = data.Poly.gpos(preLHstats.idxc);
else
  cur_pos  = preLHstats.pos_ext;
%   cur_gpos = preLHstats.gpos_ext;
end




[preLHstats.gSWj, preLHstats.gBSj, preLHstats.SWparams, preLHstats.BSparams] = composePredictedDiversity2(SWbase, SWgpos, BSbase, [], config, cur_pos, preLHstats.EgMutDiv);


bigo = 7;
