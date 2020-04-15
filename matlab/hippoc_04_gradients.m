addpath(genpath('/data/p_02323/hippoc/BrainSpace/matlab')) 
ddir      = '/data/p_02323/hippoc/data/';            % data dir
glassdir  = fullfile(ddir, 'glasserTimeseries/');    % cortex t-series
hippdir   = fullfile(ddir, 'smoothTimeseries/');     % hippocampus t-series
subjlist1 = fullfile(ddir, 'subjectListUR1QC.txt');  % 132 subjects
subjlist2 = fullfile(ddir, 'subjectListMT1QC.txt');  % 85 subjects

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

% it is HCP data, so there will be 4 scans
scans = {'rfMRI_REST1_LR', 'rfMRI_REST1_RL', ... 
    'rfMRI_REST2_LR', 'rfMRI_REST2_RL'};

% hippocampus subfields
roi_sub   = {'L_SUB', 'L_CA', 'L_DG'};
len       = [];

% here we go
for i = 1:length(ID)
    subj_glass_file = strcat(glassdir, ID{i}, '_glasserTimeseries.mat');
    subj_hipp_file  = strcat(hippdir, ID{i}, '_smoothTimeseries.mat');
    
    for j = 1:length(scans)
        % arrays
        subj_glass  = load(subj_glass_file).(scans{j});      % (1200 x 360)
        subj_hipp   = load(subj_hipp_file).(scans{j}); 
        
        % concatenate...
        for m = 1:length(roi_sub)
            subj_roi = subj_hipp.(roi_sub{m}); 
            len.(roi_sub{m}) = size(subj_roi,2);

            if m == 1
             subj_hall = subj_roi;
            else
             subj_hall = cat(2, subj_hall, subj_roi);
            end
        end

        % get correlation & sum
        A = corr(subj_hall, subj_glass);    
        if i == 1                                  % for the first subject
            SUM_CORR = A;
        else
            SUM_CORR = SUM_CORR + A;
        end
        
        fprintf('%s %s maxcorr  %.2f \n', ...
                ID{i}, scans{j}, max(max(SUM_CORR))); 
    end
end


totnum = length(ID) * length(scans)
SUM_CORR = SUM_CORR / totnum;

save('/data/p_02323/hippoc/hippocampus/matlab/avecorr_allhipsubfields.mat', 'SUM_CORR' );


% get gradients
gm = GradientMaps();
gm = gm.fit(SUM_CORR);
B = gm.gradients{1}(:,1) ; 

% assign gradient separately 
G = [];
G.L_SUB = B(1:len.L_SUB, :);
G.L_CA  = B(len.L_SUB + 1: len.L_SUB + len.L_CA );
G.L_DG  = B(len.L_SUB + len.L_CA + 1: len.L_SUB + len.L_CA + len.L_DG);


a = load('/data/p_02323/hippoc/hippocampus/matlab/surf_lsub.mat');
a = a.ave_lsub;


f = figure;
K = a;
x = (K.coord(1,:))';
y = (K.coord(2,:))';
z = -(K.coord(3,:))';
h = trisurf(K.tri,x,y,z, G.L_SUB);
set(h,'edgecolor','none')
colorbar()
caxis([-0.15, 0.15])

