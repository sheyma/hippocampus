addpath(genpath('/data/p_02323/hippoc/BrainSpace/matlab')) 

ddir      = '/data/p_02323/hippoc/data/';            % data dir
surfdir   = fullfile(ddir, 'shellsMni/');            % surface files
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

% hippocampus subfields
roi_sub   = {'L_SUB', 'L_CA', 'L_DG'};
len       = [];

% here we go

ave_lsub = [];
ave_lca  = [];
ave_dg   = [];

for i = 1:length(ID)
    
    lsub_name = strcat(surfdir, ID{i}, '_L_SUB.obj');
    lca_name  = strcat(surfdir, ID{i}, '_L_CA.obj');
    ldg_name  = strcat(surfdir, ID{i}, '_L_DG.obj');
    
    lsub = convert_surface(lsub_name);
    lca  = convert_surface(lca_name);
    ldg  = convert_surface(ldg_name);
    
    if i == 1   % first subject            
        ave_lsub.tri   = lsub.tri;
        ave_lca.tri    = lca.tri;
        ave_ldg.tri    = ldg.tri;

        ave_lsub.coord = lsub.coord;
        ave_lca.coord  = lca.coord;
        ave_ldg.coord  = ldg.coord;        
    else
        ave_lsub.tri   = ave_lsub.tri + lsub.tri;
        ave_lca.tri   = ave_lca.tri + lca.tri;
        ave_ldg.tri   = ave_ldg.tri + ldg.tri;
        
        ave_lsub.coord = ave_lsub.coord + lsub.coord;
        ave_lca.coord = ave_lca.coord + lca.coord;
        ave_ldg.coord = ave_ldg.coord + ldg.coord;
    end
   
end

totnum = length(ID);

ave_lsub.tri   = ave_lsub.tri / totnum ;
ave_lsub.coord = ave_lsub.coord / totnum;

ave_lca.tri   = ave_lca.tri / totnum ;
ave_lca.coord = ave_lca.coord / totnum;

ave_ldg.tri   = ave_ldg.tri / totnum ;
ave_ldg.coord = ave_ldg.coord / totnum;

%save('/data/p_02323/hippoc/hippocampus/matlab/surf_lsub.mat', 'ave_lsub')
%save('/data/p_02323/hippoc/hippocampus/matlab/surf_lca.mat', 'ave_lca')
%save('/data/p_02323/hippoc/hippocampus/matlab/surf_ldg.mat', 'ave_ldg')

f = figure;
K = ave_lsub;
x = (K.coord(1,:))';
y = (K.coord(2,:))';
z = -(K.coord(3,:))';

trisurf(K.tri,x,y,z,ones(1,length(x)));
colorbar()

