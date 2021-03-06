
function [thinned_pos, thinned_idx, needed] = thinnedToGrid( pos, grid )

[pos, idx] = unique(pos);

delta   = grid(1);
pos_min = grid(2);
pos_max = grid(3);


% needed = zeros(size(pos));
% 
% for curpos = pos_min:delta:pos_max
%   needed(find(pos>=curpos,1,'first')) = 1;
%   needed(find(pos<=curpos,1,'last' )) = 1;
% end
% thinned_pos = pos(find(needed));


needed = histc([pos_min:delta:pos_max], [-Inf; pos ;Inf]);
thinned_i1 = find(needed(1:end-2)>0);
thinned_i2 = [thinned_i1(2:end)-1  length(pos)];

thinned_pos = pos([thinned_i1 thinned_i2]);
thinned_idx = idx([thinned_i1 thinned_i2]);

[thinned_pos, idx0] = unique(thinned_pos);
thinned_idx = thinned_idx(idx0);


bigo = 7;
