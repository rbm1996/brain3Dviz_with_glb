clear all

% add functions to path
subdir = genpath("myUtils");
addpath(subdir);


%% load brain nodes data

shaeferRight = [101 : 200 208:214];


corticalIdx = shaeferRight;

infotable = readtable("Schaefer200_allinfo.csv");

xyzSub = table2array(infotable(:,7:9));

xyz = xyzSub(corticalIdx , :);

networksAssignment = table2array(infotable(:,5));
nNetworks = max(networksAssignment);

% schaefer colors
load('mySchaeferColorMap_9nets.mat')
color = myColorMap(networksAssignment(corticalIdx) , :);


%% load connectome 
subName = "2134464";
data_dir = 'myConnectomes';
tablePath = fullfile(data_dir , sprintf('%s',subName));

tableName = fullfile(tablePath,sprintf(...
    '%s_connectome_edgelist_schaefer200-yeo17_rmap_s1dild.csv' , ...
    subName )) ;
A = load_UKBiobank_connectome(tableName , 216);


A = A/max(max(A));
allStr = sum(A);
strengthNode = 2* sqrt(sqrt(allStr(corticalIdx)));




%% generate spheres

strengthNode = 1.6 * strengthNode.^0.2;

alphaBrain = 0.8;

radius_factor = 0.5;
distance_factor = 1;

nodes = generate_nodes(strengthNode, xyz, radius_factor, color, distance_factor);

nodes.vertices = nodes.vertices - ones(size(nodes.vertices , 1) ,1)* mean(nodes.vertices);
nodes.vertices(:,2) = nodes.vertices(:,2) + 5;
nodes.vertices(:,1) = nodes.vertices(:,1) - 2;


% center
% brain.vertices = brain.vertices - ones(size(brain.vertices , 1) ,1)* mean(brain.vertices);
%
% tmp_shading_color = brain.shading_pre * diff(range);
% tmp_shading_color = tmp_shading_color - min(tmp_shading_color) + range(1);
% shading_color = repmat(tmp_shading_color,1,3);
% brain.color = shading_color;

% Parameters
range = [0 1];

nodes.colors = nodes.colors * diff(range);
nodes.xyz = nodes.xyz - ones(size(nodes.xyz , 1) ,1)* mean(nodes.vertices);


%% correct offset xyz


nodes.vertices = nodes.vertices - ones(size(nodes.vertices , 1) ,1)* mean(nodes.vertices);

coordToShift = 1; offSet = -30;
nodes.vertices(:,coordToShift) = nodes.vertices(:,coordToShift) - offSet;

coordToShift = 2; offSetY = 0.5;
nodes.vertices(:,coordToShift) = nodes.vertices(:,coordToShift) - offSetY;


xyz = xyz - ones(size(xyz , 1) ,1)* mean(xyz);

coordToShift = 1; offSet = -30;
xyz(:,coordToShift) = xyz(:,coordToShift) - offSet;

coordToShift = 2; offSetY = 0.5;
xyz(:,coordToShift) = xyz(:,coordToShift) - offSetY;



%% Plot left hemisphere surface and right hemisphere spheres

figure;


alphaBrain = 1;
range = [0 0.7];

h_brain = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaBrain , 'lh' , myColorMap(networksAssignment,:)) ;


hold on
spheresPlot = trisurf(nodes.indices , ...
    nodes.vertices(:,1),nodes.vertices(:,2),nodes.vertices(:,3), ...
    'edgecolor','none', 'FaceLighting', 'gouraud', 'AmbientStrength', 0.5, ...
    'FaceVertexCData', nodes.colors);


% trisurf(allIndices, allVertices(:, 1), allVertices(:, 2), allVertices(:, 3), C, 'edgecolor','none');
% light
% shading interp
shading interp
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
grid off
axis off
view(0, 90);

%%  get links from connectome
nNodes = length(corticalIdx);

adj = A(corticalIdx , corticalIdx);

adj = adj - diag(diag(adj));
% adj = triu(adj);

% th = min(min(adj)) + (1 - 0.85) * (max(max(adj)) - min(min(adj)));
% adj(adj < th) = 0;

% adj(101:107 , 101:107) = 0;

% manual adjustment to remove some links
adj(101:107 , :) = 0;
adj(: , 101:107) = 0;

adj(13 , 101) = 0.2;
adj(101 , 102) = 1;
adj(103 , 102) = 1;
adj(101 , 103) = 1;

adj(103 , 104) = 1;
adj(104 , 105) = 1;
adj(105 , 106) = 1;
adj(106 , 107) = 0.2;

adj(105 , 101) = 1;

adj(101:107 , 101:107) = (adj(101:107 , 101:107) + adj(101:107 , 101:107)');


%%% extract backbone of the network
[CIJtree,CIJclus] = backbone_wu(adj,6);
CIJtreeTri = triu(CIJtree);

% find links
[row , col , wei] = find(CIJtreeTri);

myLinks = [row col];


% myEdgesWidth = 0.5 * ones(1,size(myLinks , 1)) ;


linkMetric = (wei).^0.3;

linkMetric = (linkMetric - min(linkMetric))/ (max(linkMetric) - min(linkMetric)) ;

myEdgesWidth = 0.1 + 2 * linkMetric;

links = generate_links_v2(strengthNode, myLinks , myEdgesWidth , xyz , radius_factor, color, distance_factor);

% links.vertices = links.vertices - ones(size(links.vertices , 1) ,1)* mean(links.vertices);

%% plot the links
hold on


patch('Faces',links.faces,'Vertices',links.vertices, 'FaceVertexCData' ,links.colors , 'FaceColor','interp' , 'EdgeColor' , 'none');



%% stack all patches : stack vertices, faces , colors
ax = gca;ch = ax.Children; % get axes and children 
nPatches = length(ch);

nFaces = 0;
nVert = 0;

allFaces = [];
allVert = [];
allVertColors = [];

kPatch = 1;
allFaces = [allFaces ; ch(kPatch).Faces(:, [1 2 3]) + nVert];
allFaces = [allFaces ; ch(kPatch).Faces(:, [1 3 4]) + nVert];

allVert = [allVert ; ch(kPatch).Vertices];
allVertColors= [allVertColors ; ch(kPatch).FaceVertexCData];
nVert = nVert + size(ch(kPatch).Vertices , 1);

for kPatch = [nPatches-1 nPatches]
    allFaces = [allFaces ; ch(kPatch).Faces + nVert];
    allVert = [allVert ; ch(kPatch).Vertices];
    allVertColors= [allVertColors ; ch(kPatch).FaceVertexCData];
    nVert = nVert + size(ch(kPatch).Vertices , 1);

end


%% save GLB object
% https://github.com/dmitrishastin/matlab2glb
% look the website for details on the parametrs
% generate .glb

example.POSITION = allVert; % vertices
example.indices = allFaces;% faces
example.COLOR_0 = allVertColors.^2.2; % correct colors for glb format scale
aa = 0.55;
example.prop.material.pbrMetallicRoughness.baseColorFactor =  [aa aa aa 1];
example.prop.material.pbrMetallicRoughness.metallicFactor = 0.5; % 1
example.prop.material.pbrMetallicRoughness.roughnessFactor = 1; % 0.1
example.prop.material.alphaMode = 'OPAQUE' ; % 'BLEND';
example.prop.material.doubleSided = true;

% save file
write_glb('myBrain_3D.glb', example);

disp('Done!')


