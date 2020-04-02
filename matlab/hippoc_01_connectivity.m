% definde data directories

ddir     = '/data/p_02323/hippoc/data/';            % data dir

glassdir = fullfile(ddir, 'glasserTimeseries/');    % cortex t-series
hippdir  = fullfile(ddir, 'smoothTimeseries/');     % hippocampus t-series
subjlist = fullfile(ddir, 'subjectListUR1QC.txt');  % 132 subjects

% it is HCP data, so there will be 4 scans

scans = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL', ... 
    'rfMRI_REST2_LR', 'rfMRI_REST2_RL'};

% hippocampus segmentations

roi = {'L_SUB', 'R_SUB', ...
       'L_CA', 'R_CA', ...
       'L_DG', 'R_DG'};

% get subject ID as type 'cell'

fid      = fopen(subjlist); 
C        = textscan(fid,'%s', 'CollectOutput',1);  
fclose(fid);
ID       = C{1}(:,1); 

% get subibuculum(hippocampus)-to-cortex(glasser) connectivity

roi_sub = {'L_SUB', 'R_SUB'};

glasservertexnum = 360; 
SUM_CORR         = zeros(glasservertexnum, 1); 
totnum           = 0;

for i = 1:length(ID)
    for j = 1:length(scans)
        for m = 1:length(roi_sub)

            % files
            subj_glass_file = strcat(glassdir, ID{i}, '_glasserTimeseries.mat');
            subj_hipp_file  = strcat(hippdir, ID{i}, '_smoothTimeseries.mat');

            % arrays
            subj_glass  = load(subj_glass_file).(scans{j});              % (1200 x 360)
            subj_hipp   = load(subj_hipp_file).(scans{j}).(roi_sub{m});  % (1200 x 1024)
            subj_hippav = mean(subj_hipp, 2);                            % (1200 x 1)
            subj_corr   = corr(subj_glass, subj_hippav);                 % (360 x 1)

            fprintf('%s %s %s maxcorr  %.2f \n', ...
                    ID{i},  roi_sub{m}, scans{j}, max(subj_corr));

            % Fisher RtoZ prior to averaging
            SUM_CORR = SUM_CORR + atanh(subj_corr); 
            totnum   = totnum + 1;

        end
    end
end

SUM_CORR = SUM_CORR / totnum;

% plotting... (install github repo's BrainSpace and gifti-master)

addpath(genpath('/data/p_02323/hippoc/BrainSpace/matlab')) % plotting tool
addpath(genpath('/data/p_02323/hippoc/gifti-master/'))     % gifti tool

labeling_glasser= load(fullfile(ddir, 'glasser.csv'));     % 64k labeling

[surf_lh, surf_rh] = load_conte69();      % 32k left & 32k right fsaverage

% here we go
obj = plot_hemispheres(SUM_CORR, {surf_lh,surf_rh}, ...
                       'views','lmap', ...
                       'labeltext', {'<SUBICULUM>'}, ...
                       'parcellation', labeling_glasser);

obj.colorlimits([0, 0.5]);
obj.labels('FontSize',15)
colormap(obj.handles.figure, hot)

