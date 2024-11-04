function links = generate_links_v2(Xradius, myLinks , myEdgesWidth , xyz, radius_factor, color, distance_factor, varargin)
% nodes positions
% addpath('/Users/juliana.gonzalez/ownCloud/graph_analysis/');
% results_file = 'resultsROI_Subject001_Condition001.mat';
% results = load(results_file);
% xyz = results.xyz;

%addpath('/Users/juliana.gonzalez/ownCloud/github/synesnet/');
%X_file = 'Xnet_syn_mean.mat';
%load(X_file);


%% new for links
numLinks = size(myLinks,1);
% Initialize arrays for storing vertices and indices
links.vertices = [];
links.indices = [];
links.colors = [];


% %% rotation matrices
% syms t
% 
% Rx = [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)];
% Ry = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
% Rz = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];
% 

R = 0.5;
Nr=10;
Nt=10;
zz=0:0.2:1;

links = struct;
links.vertices = [];
links.faces = [];
links.colors = [];


%% loop over links
nVerticesTotal = 0;

myGray  = gray();

offsetColor = 70;
myGray = [ones(offsetColor , 1)*myGray(1,:) ; myGray(1 : end - offsetColor,:)];

flipGray = flip(myGray);

for kLink = 1:numLinks
    fprintf('link %i / %i \n' , kLink , numLinks)
    node1 = xyz(myLinks(kLink , 1) , :);
    node2 = xyz( myLinks(kLink , 2) , :);

    %% generate original cylinder


    Nr = 50;
    Nt =20;

    R = myEdgesWidth(kLink)* ones(1,Nt) / 2;

    [X, Y, Z] = cylinder2P(R, Nr,node1,node2);

    % figure;
    % surObj = surf(X,Y,Z);

    fvc = surf2patch(surf(X,Y,Z));

    cData = fvc.facevertexcdata;


    %%   %%%%%
    
    nVertices = size(fvc.vertices , 1);

    links.vertices = [links.vertices ; fvc.vertices];
    links.faces = [links.faces ; fvc.faces + nVerticesTotal];
    links.colors = [links.colors ; flipGray(1 + floor(255*(cData-min(cData))/(max(cData)-min(cData))) , :)];

    nVerticesTotal = nVerticesTotal+ nVertices; 

end