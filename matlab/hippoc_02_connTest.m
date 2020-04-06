
%%%%% set fsaverage coordinates

addpath(genpath('/data/p_02323/hippoc/BrainSpace/matlab')) % plotting tool
addpath(genpath('/data/p_02323/hippoc/gifti-master/'))     % gifti tool
addpath(genpath('/data/p_02323/hippoc/micaopen'))


[surf_lh, surf_rh] = load_conte69();      % 32k left & 32k right fsaverage

D = [];
D.coord = [surf_lh.coord, surf_rh.coord];
D.tri   = [surf_lh.tri; surf_rh.tri + length(surf_lh.coord)];

length(D.coord)                           % 64k

BoSurfStatViewData(rand(length(D.coord),1), D, 'random brain')
colormap('hot')

%%%%% get subject-specific connectivity
ddir     = '/data/p_02323/hippoc/data/';           
glassdir = fullfile(ddir, 'glasserTimeseries/');    % cortex t-series
hippdir  = fullfile(ddir, 'smoothTimeseries/');     % hippocampus t-series
subjlist = fullfile(ddir, 'subjectListUR1QC.txt');  % 132 subjects

scans   = {'rfMRI_REST1_LR'};
roi_sub = {'L_SUB', 'R_SUB'};

fid      = fopen(subjlist); 
txt      = textscan(fid,'%s', 'CollectOutput',1);  
fclose(fid);
ID       = txt{1}(:,1); 

C360_all = zeros(length(ID), 360);

for i = 1:length(ID)
    
    k = zeros(360, 1);
    
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
            
            k = k + atanh(subj_corr);
        
        end
    end
    
    k = k / (length(scans) * length(roi_sub)); 
    fprintf('%s  maxcorr  %.2f \n', ID{i}, max(k));    
    
    C360_all(i,:) = k; 
end


%%%%% resample 360 -->> 64k for each subject & plot mean connectivity

mylabel  = load(fullfile(ddir, 'glasser.csv'));     % 64k labeling

C64k_all = zeros(length(ID), 64984);                

for i = 1:length(ID)
    for j = 1:360
       C64k_all(i, (find(mylabel == j))) = C360_all(i, j); 
    end
end

BoSurfStatViewData(mean(C64k_all, 1), D, 'average connectivity')
BoSurfStatColLim([0 0.5])
colormap('hot')

%%%%% one-sample t-test across subjects

T        = C64k_all;                    
subjID   = ID;
contrast = ones(length(ID),1);
M        = 1 + term(contrast); 
slm      = SurfStatLinModS(T, M, D); 
slm      = SurfStatT(slm, contrast);

Tvals    = slm.t;
Tvals(Tvals < 15) = Inf;                          % thresholding
BoSurfStatViewData(Tvals, D, 't-values')
BoSurfStatColLim([15 40])
colormap([hot; .7 .7 .7])

