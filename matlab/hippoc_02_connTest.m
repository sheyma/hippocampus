
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
ddir      = '/data/p_02323/hippoc/data/';           
glassdir  = fullfile(ddir, 'glasserTimeseries/');    % cortex t-series
hippdir   = fullfile(ddir, 'smoothTimeseries/');     % hippocampus t-series
subjlist1 = fullfile(ddir, 'subjectListUR1QC.txt');  % 132 subjects
subjlist2 = fullfile(ddir, 'subjectListMT1QC.txt');  % 85 subjects

scans = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL', ... 
    'rfMRI_REST2_LR', 'rfMRI_REST2_RL'};

roi_sub = {'L_SUB', 'R_SUB'};

% get subject ID's as cell

fid      = fopen(subjlist1); 
txt      = textscan(fid,'%s', 'CollectOutput',1);  
fclose(fid);
ID1      = txt{1}(:,1); 

fid      = fopen(subjlist2); 
txt      = textscan(fid,'%s', 'CollectOutput',1);  
fclose(fid);
ID2      = txt{1}(:,1); 

ID = [ID1; ID2];

% here we go...

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

C360_all = load('/data/p_02323/hippoc/hippocampus/matlab/avecorr_217.mat');
C360_all = C360_all.C360_all; 

%%%% plot average connectivity across connectivity

C360_ave  = mean(C360_all,1);
C360_surf = zeros(64984, 1);   

mylabel   = load(fullfile(ddir, 'glasser.csv'));     % 64k labeling

for i = 1:360    
    C360_surf(mylabel == i) = C360_ave(i); 
end

BoSurfStatViewData(C360_surf, D, 'average connectivity')
BoSurfStatColLim([0 0.5])
colormap('hot')

%%%% one-sample t-test

parcels  = C360_all;                % 217 x 360
subjID   = ID;                      % 217 
contrast = ones(length(subjID),1);  % 217 x 1
M        = 1 + term(contrast);      % 217 x 1

slm = SurfStatLinMod(parcels, M);
slm = SurfStatT(slm, contrast);

Tvals = slm.t;

% multiple comparison correction: Benferroni
pvals = 1-tcdf(slm.t, slm.df);
pvals = pvals*size(pvals,2);

Tsurf = zeros(64984, 1);   
Psurf = zeros(64984, 1);  

for i = 1:360    
    Tsurf(mylabel == i) = Tvals(i); 
    Psurf(mylabel == i) = pvals(i); 

end

Tsurf(Tsurf < 20) = Inf;                          % thresholding
BoSurfStatViewData(Tsurf, D, 't-values')
BoSurfStatColLim([20 70])
colormap([hot; .7 .7 .7])
                   
BoSurfStatViewData(Psurf, D, 'p-values')
BoSurfStatColLim([0 0.05])
colormap([parula; .7 .7 .7])

