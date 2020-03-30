

ddir     = '/data/p_02323/hippoc/data/';            % data dir

glassdir = fullfile(ddir, 'glasserTimeseries/');    % cortex t-series
hippdir  = fullfile(ddir, 'smoothTimeseries/');     % hippocampus t-series
subjlist = fullfile(ddir, 'subjectListUR1QC.txt');  % 132 subjects


scans = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL', ... 
    'rfMRI_REST2_LR', 'rfMRI_REST2_RL'};

roi = {'L_SUB', 'R_SUB',
    'L_CA', 'R_CA',
    'L_DG', 'R_DG'};

% get subject ID as type 'cell'

fid      = fopen(subjlist); 
C        = textscan(fid,'%s', 'CollectOutput',1);  
fclose(fid);
ID       = C{1}(:,1); 


% example correlation

i = 1;          % ID{1}    : 'HCP_100610'
j = 1;          % scans{1} : 'rfMRI_REST1_LR'
m = 1;          % roi{1}   : 'L_SUB'


subj_glass_file = strcat(glassdir, ID{i}, '_glasserTimeseries.mat')
subj_hipp_file  = strcat(hippdir, ID{i}, '_smoothTimeseries.mat')

subj_glass  = load(subj_glass_file).(scans{j});         % (1200 x 360)

subj_hipp   = load(subj_hipp_file).(scans{j}).(roi{m}); % (1200 x 1024)
subj_hippav = mean(subj_hipp, 2);                       % (1200 x 1)


a = corr(subj_glass, subj_hippav);                      % (360 x 1)

