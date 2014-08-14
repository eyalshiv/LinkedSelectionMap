
function fdata = LS_LoadVariationData( inputfiles ) %inputfiles_file

if isstr(inputfiles)
  inputfiles = file2struct( inputfiles );
end
  

%% general inits

% [ff_chr_id, ff_chr_len] = textread(chr_features_file, '%s\t%d', 'headerlines', 1);

C = length(inputfiles.poly);

for c=1:C
  f = fopen( inputfiles.poly{c}, 'rt' );
  line = fgetl(f);
  [fdata{c}.chr_id, fdata{c}.chr_len] = sscanf( line, '#%s\t%d' );
  f = fclose(f);
end

% for c=1:C
%   ichr(c) = find(strcmp( ff_chr_id, chr_id{c} ));
%   chr_len(c) = ff_chr_len(ichr(c));
% end


%% load genetic map

for c=1:C
  [genmap{c}.pos, genmap{c}.c, genmap{c}.R]  = textread( inputfiles.genmap{c}, '%d\t%f\t%f', 'headerlines', 1 );
  genmap{c}.file = inputfiles.genmap{c};
  genmap{c}.name = inputfiles.genmap_token;
  
  fdata{c}.genmap = genmap{c};
  fdata{c}.genmaplims = [min(genmap{c}.R) max(genmap{c}.R)] /100;
end



%% load poly data

for c=1:C
  f = fopen( inputfiles.poly{c}, 'rt' );

% [pos strand sample freq anc der polyType]
  Z = textscan( f, '%d\t%c\t%d\t%f\t%s\t%s\t%d', 'commentstyle', '#' );

  L = length(Z{1});
  
  fdata{c}.Poly.pos  = Z{1};
  fdata{c}.Poly.gpos = applyGenmap2pos( fdata{c}.genmap, fdata{c}.Poly.pos );
  fdata{c}.Poly.strand  = Z{2};
  
  idxS    = find( Z{7} == 1 ); % 0=mono, 1=SYN, 2=NS
  fdata{c}.Poly.idxS    = idxS;
  fdata{c}.Poly.sample  = Z{3};
  fdata{c}.Poly.specS   = [Z{4}(idxS) double(Z{3}(idxS))-Z{4}(idxS)];
  fdata{c}.Poly.anc = char(zeros([L 3]));
  fdata{c}.Poly.der = char(zeros([L 3]));
  for i=1:L
  fdata{c}.Poly.anc(i,:)     = Z{5}{i};
  fdata{c}.Poly.der(i,:)     = Z{6}{i};
  end
%   fdata{c}.Poly.type    = Z{7};
  
  f = fclose(f);
end



%% load mutation rate proxy data

for c=1:C
  f = fopen( inputfiles.mutprox{c}, 'rt' );

% [pos support mutrate]
  Z = textscan( f, '%d\t%f\t%f', 'commentstyle', '#' );

  
  fdata{c}.MutProx.pos  = Z{1};
  fdata{c}.MutProx.gpos = applyGenmap2pos( fdata{c}.genmap, fdata{c}.MutProx.pos );
  
  fdata{c}.MutProx.paml_codons  = Z{2}; 
  fdata{c}.MutProx.dS       = Z{3};
  fdata{c}.MutProx.dS( fdata{c}.MutProx.paml_codons<=0 ) = -1;
  
  f = fclose(f);
end

  
end
