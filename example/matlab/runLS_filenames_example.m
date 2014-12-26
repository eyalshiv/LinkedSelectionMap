


% root dir
base_dir = 'E:/GitHub/LinkedSelectionMap/example';




%% chromosomes

% a file containing a list of chromosomes and the length of each one in bp
chr_features_file = [base_dir '/data/chromosomal_features_example.txt'];
[chr_id, chr_len] = textread(chr_features_file, '%s %d', 'headerlines', 1);

% the number of chromosomes
C = length(chr_id);


%% genetic maps

genmap_token     = 'Comeron';
for c=1:C
  % genetic map file names for each chromosome
  genmap_files{c}     =  [base_dir sprintf('/data/genetic_maps/genmap_Comeron_%s.txt',chr_id{c})];
end



%% selection annotations and LS maps

% enumerate CS annotations
iSWAnnoExonicNS = 1;
iSWAnnoUTR = 2;
iSWAnnoIntronic = 3;

% label CS annotations
SW_anno_tokens    = {...
  'exonicNS',...
  'UTR',...
  'intronic',...
  '' };
%   'intergenic' };

% prefices of CS annotation (i.e. substitutions list) files
SW_anno_fileprefs = {...
  [base_dir '/data/annotations/substitutions_exonicNS_'],...
  [base_dir '/data/annotations/substitutions_UTR_'],...
  [base_dir '/data/annotations/substitutions_intronic_'],...
  '' };
% we have no data for substitutions at intergenic regions
%  [base_dir '/data/annotations/substitutions_intergenic_'] };

% enumerate BS annotations
iBSAnnoExonic = 1;
iBSAnnoUTR = 2;
iBSAnnoIntronic = 3;
iBSAnnoIntergenic = 4;

% label BS annotations
BS_anno_tokens    = {...
  'exons',...
  'UTRs',...
  'longintrons',...
  'intergenic'};

% prefices of BS annotation (i.e. conserved segments) files
BS_anno_fileprefs = {...
  [base_dir '/data/annotations/longest_transcript_533_exon_'],...
  [base_dir '/data/annotations/longest_transcript_533_UTR_'],...
  [base_dir '/data/annotations/longest_transcript_533_longintron_'],...
  [base_dir '/data/annotations/longest_transcript_533_intergenic_'] };

% prepare the file names of all annotations across all chromosomes
for c=1:C
%   for b=1:length(BS_anno_tokens)
    BS_anno1_files{c} = '';
    BS_anno2_files{c} = '';
    BS_anno3_files{c} = '';
    BS_anno4_files{c} = '';
    if ~isempty(BS_anno_fileprefs{1}), BS_anno1_files{c} = [BS_anno_fileprefs{1} chr_id{c} '.coords'];, end
    if ~isempty(BS_anno_fileprefs{2}), BS_anno2_files{c} = [BS_anno_fileprefs{2} chr_id{c} '.coords'];, end
    if ~isempty(BS_anno_fileprefs{3}), BS_anno3_files{c} = [BS_anno_fileprefs{3} chr_id{c} '.coords'];, end
    if ~isempty(BS_anno_fileprefs{4}), BS_anno4_files{c} = [BS_anno_fileprefs{4} chr_id{c} '.coords'];, end
%   end
%   for b=1:length(SW_anno_tokens)
    SW_anno1_files{c} = '';
    SW_anno2_files{c} = '';
    SW_anno3_files{c} = '';
    SW_anno4_files{c} = '';
    if ~isempty(SW_anno_fileprefs{1}), SW_anno1_files{c} = [base_dir sprintf('/data/annotations/substitutions_%s_%s.txt', chr_id{c}, SW_anno_tokens{1})];, end
    if ~isempty(SW_anno_fileprefs{2}), SW_anno2_files{c} = [base_dir sprintf('/data/annotations/substitutions_%s_%s.txt', chr_id{c}, SW_anno_tokens{2})];, end
    if ~isempty(SW_anno_fileprefs{3}), SW_anno3_files{c} = [base_dir sprintf('/data/annotations/substitutions_%s_%s.txt', chr_id{c}, SW_anno_tokens{3})];, end
    if ~isempty(SW_anno_fileprefs{4}), SW_anno4_files{c} = [base_dir sprintf('/data/annotations/substitutions_%s_%s.txt', chr_id{c}, SW_anno_tokens{4})];, end
    %     SW_anno_files{c}{b} = [SW_anno_fileprefs{b} chr_id{c} '.txt'];
%   end
end

% prepare CS grid file names (these are lists of positions at which to calculate the efects of CS, instead of using a fixed-distance grid)
for c=1:C
  SWbase_pos_grid_files{c} = [base_dir '/data/annotations/SWbase_positions_' chr_id{c} '.txt'];
end



%% variation data

% for c=1:C
%   poly_files{c} = { [base_dir sprintf('/data/Bailor/strains/Bailor_%s_output1.txt', chr_id{c})],...
%                     [base_dir sprintf('/data/Bailor/strains/Bailor_%s_output2.txt', chr_id{c})],...
%                     [base_dir sprintf('/data/Bailor/strains/Bailor_%s_output3.txt', chr_id{c})]};
% end

% prepare file names of polymorphism data and approximated local mutation
% rates
for c=1:C
  poly_proc_files{c}       = [base_dir sprintf('/data/neutral_variation/Filtered_140521_poly_%s.txt',    chr_id{c})];
  MutProx_proc_files{c}    = [base_dir sprintf('/data/neutral_variation/Filtered_140521_mutrate_%s.txt', chr_id{c})];
end

% prepare mask files for the inference and for goodness-of-fit evaluation
for c=1:C
  inference_mask_files{c}       = [base_dir sprintf('/work/masks/codons_mask_inference_%s.txt',  chr_id{c})];
  evaluation_mask_files{c}      = [base_dir sprintf('/work/masks/codons_mask_evaluation_%s.txt', chr_id{c})];
%   bootstrap_output_mask_pref{c} = [base_dir sprintf('/work/bootstrap/masks/codons_mask_inference_%s', chr_id{c})];
end


% wrap all relevant input and output filenames in a configuration struct
infcfg_file = [base_dir '/work/example_cfg.txt'];

files_invar.poly           = poly_proc_files;
files_invar.mutprox        = MutProx_proc_files;
files_invar.genmap         = genmap_files;
files_invar.genmap_token   = genmap_token;
files_invar_file           = [base_dir '/work/files_invar.txt'];

files_buildGE.outdir            = [base_dir '/work/GEs/']; 
files_buildGE.outfile_pref      = '';
files_buildGE.chr_features_file = chr_features_file;
files_buildGE.chr_id            = chr_id; %(1:4)
files_buildGE.pos_grid_files    = SWbase_pos_grid_files; %(1:4)
files_buildGE.genmap_files      = genmap_files; %(1:4)
files_buildGE.genmap_token      = genmap_token;
files_buildGE.SW_anno1_files = SW_anno1_files; %(1:4)
files_buildGE.SW_anno2_files = SW_anno2_files; %(1:4)
files_buildGE.SW_anno3_files = SW_anno3_files; %(1:4)
files_buildGE.SW_anno4_files = SW_anno4_files; %(1:4)
files_buildGE.SW_anno_tokens = SW_anno_tokens;
files_buildGE.BS_anno1_files = BS_anno1_files; %(1:4)
files_buildGE.BS_anno2_files = BS_anno2_files; %(1:4)
files_buildGE.BS_anno3_files = BS_anno3_files; %(1:4)
files_buildGE.BS_anno4_files = BS_anno4_files; %(1:4)
files_buildGE.BS_anno_tokens = BS_anno_tokens;
files_buildGE_file           = [base_dir '/work/files_buildGE.txt'];


files_masks.inference           = inference_mask_files;
files_masks.evaluation          = evaluation_mask_files;
files_masks_file                = [base_dir '/work/files_masks.txt'];


files_bootstrap.chr_features_file = chr_features_file;
files_bootstrap.input_mask        = files_masks.inference;
% files_bootstrap.output_mask_pref  = bootstrap_output_mask_pref;
% files_bootstrap_file              = [base_dir '/work/bootstrap/files_bootstrap_masks.txt'];



