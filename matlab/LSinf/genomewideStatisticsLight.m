
function gwStats = genomewideStatisticsLight( fdata )


if nargin < 2
  poolNcorrStats = 0;
end


Ls2L = 0.25;
Ln2L = 1-Ls2L;



C = length(fdata);

codonsSyn1 = 0; codonsDat1 = 0; codonsSyn2 = 0; codonsDat2 = 0; codonsSyn3 = 0; codonsDat3 = 0; codonsSyn3b = 0; codonsDat3b = 0; codonsSyn4 = 0; codonsDat4 = 0;  codonsSyn5 = 0; codonsDat5 = 0;
gMutDiv = [];

% gwStats.nNSsubs         = 0;
% gwStats.nNSsubsfake     = 0;
gwStats.mapL            = 0;

for c=1:C
  %     codonsSyn1 = codonsSyn1 + length(fdata{c}.LongDiv.idxS);
  %     codonsDat1 = codonsDat1 + length(fdata{c}.LongDiv.gpos);
  %     codonsSyn2 = codonsSyn2 + length(fdata{c}.Div.idxS);
  %     codonsDat2 = codonsDat2 + length(fdata{c}.Div.gpos(fdata{c}.Div.isdat));
  codonsSyn3 = codonsSyn3 + length(fdata{c}.Poly.idxS);
  codonsSyn3b= codonsSyn3b+ sum(prod(fdata{c}.Poly.specS,2));
  codonsDat3 = codonsDat3 + length(fdata{c}.Poly.gpos);
  codonsDat3b= codonsDat3b+ sum(fdata{c}.Poly.sample.*(fdata{c}.Poly.sample-1)/2);
  
  %     codonsSyn4 = codonsSyn4 + length(intersect(fdata{c}.LongDiv.gpos(fdata{c}.LongDiv.idxS), fdata{c}.Poly.gpos));
  %     codonsDat4 = codonsDat4 + length(intersect(fdata{c}.LongDiv.gpos,                        fdata{c}.Poly.gpos));
  %     codonsSyn5 = codonsSyn5 + length(intersect(fdata{c}.Div.gpos(fdata{c}.Div.idxS),  fdata{c}.Poly.gpos));
  %     codonsDat5 = codonsDat5 + length(intersect(fdata{c}.Div.gpos(fdata{c}.Div.isdat), fdata{c}.Poly.gpos));
  
  pos = double(fdata{c}.Poly.pos);
  dSidx = find( ~isnan(fdata{c}.MutProx.dS) | fdata{c}.MutProx.dS < 0 ); 
  cgMutDiv   = interp1(double(fdata{c}.MutProx.pos(dSidx)), fdata{c}.MutProx.dS(dSidx), pos);
  cgMutDiv(pos < fdata{c}.MutProx.pos(dSidx(1  ))) = fdata{c}.MutProx.dS(dSidx(1  ));
  cgMutDiv(pos > fdata{c}.MutProx.pos(dSidx(end))) = fdata{c}.MutProx.dS(dSidx(end));
  gMutDiv = [gMutDiv; cgMutDiv];
  
  gwStats.mapL = gwStats.mapL + max(fdata{c}.Poly.gpos) - min(fdata{c}.Poly.gpos);
  %     gwStats.nNSsubs        = gwStats.nNSsubs     + length(fdata{c}.Div.idxN);
  %     gwStats.nNSsubsfake    = gwStats.nNSsubsfake + length(fdata{c}.Div.idxN) + length(fdata{c}.Div.idxNfake);
  
  gwStats.chr_len(c)       = fdata{c}.chr_len;
end

% gwStats.remoteDiv       = (codonsSyn1/codonsDat1);% / (3*Ls2L);
% gwStats.remoteDivAtPoly = (codonsSyn4/codonsDat4);% / (3*Ls2L);
% gwStats.closeDiv        = (codonsSyn2/codonsDat2);% / (3*Ls2L);
% gwStats.closeDivAtPoly  = (codonsSyn5/codonsDat5);% / (3*Ls2L);
gwStats.Poly            = (codonsSyn3/codonsDat3);% / (3*Ls2L);
gwStats.Het             = (codonsSyn3b/codonsDat3b);% / (3*Ls2L);
gwStats.MutProx         = nanmean(gMutDiv);% / (3*Ls2L);
  
% gwStats.remote2closeDiv = gwStats.remoteDiv / gwStats.closeDiv;



end

