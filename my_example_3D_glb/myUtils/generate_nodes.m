function nodes = generate_nodes(Xradius, xyz, radius_factor, color, distance_factor, varargin)
    % nodes positions
    % addpath('/Users/juliana.gonzalez/ownCloud/graph_analysis/');
    % results_file = 'resultsROI_Subject001_Condition001.mat';
    % results = load(results_file);
    % xyz = results.xyz;
    
    %addpath('/Users/juliana.gonzalez/ownCloud/github/synesnet/');
    %X_file = 'Xnet_syn_mean.mat';
    %load(X_file);
    
    % Initialize an empty vector for storing x positions
    xVector = [];
    yVector = [];
    zVector = [];
    
    % Iterate over each cell in the xyzCell array
    % distance_factor = 1.15;
    for i = 1:size(xyz,1)
        % Access the x position from the current cell and append it to the xVector
        xVector = [xVector, xyz(i, 1) * distance_factor] ;
        yVector = [yVector, xyz(i, 2) * distance_factor] ;
        zVector = [zVector, xyz(i, 3) * distance_factor];
    end
    
    % Parameters
    % radius_factor = 0.1;
    numSpheres = size(xyz,1)-1;       % Number of spheres
    radius = Xradius * radius_factor;  % Radii of the spheres
    
    numSegments = 20;    % Number of vertical segments
    numSlices = 40;      % Number of horizontal slices
    
    % Initialize arrays for storing vertices and indices
    nodes.vertices = [];
    nodes.indices = [];
    nodes.colors = [];
    nodes.xyz = [];
    
    % Min-max normalization of sphere sizes
    minRadius = min(radius);
    maxRadius = max(radius);
    normalizedRadius = (radius - minRadius) / (maxRadius - minRadius) + 0.1;
    normalizedRadius = normalizedRadius * 4;

    % Generate spheres
    for i = 1:numSpheres
        % Generate sphere vertices
        theta = linspace(0, 2*pi, numSlices+1);
        phi = linspace(-pi/2, pi/2, numSegments+1);
        [theta, phi] = meshgrid(theta, phi);
        x = xVector(i) + normalizedRadius(i) * cos(theta) .* cos(phi);
        y = yVector(i) + normalizedRadius(i) * sin(theta) .* cos(phi);
        z = zVector(i) + normalizedRadius(i) * sin(phi);
        vertices = [x(:), y(:), z(:)];
        nodes.vertices = [nodes.vertices; vertices];
        node_xyz = [xVector(i), yVector(i), zVector(i)];
        nodes.xyz = [nodes.xyz; node_xyz];



        % color
        n_color = color(i,1:end);
        shading = (vertices(:, 3) + normalizedRadius(i) - zVector(i)) / (2 * normalizedRadius(i));  % Vertical position-based shading
        shading = 0.5 + shading / 2;  % Lighter shadow effect (adjust the scaling factor for desired brightness)

        numVertices = size(vertices, 1);
        color_v = repmat(n_color, numVertices, 1);  
        color_v = repmat(shading, 1, 3) .* color_v;  % Apply shading to the color vector
        nodes.colors = [nodes.colors; color_v];
    
        % Generate sphere indices
        indices = [];
        for j = 1:numSegments
            for k = 1:numSlices
                k1 = (j - 1) * (numSlices + 1) + k;
                k2 = k1 + 1;
                k3 = k1 + numSlices + 2;
                k4 = k1 + numSlices + 1;
                indices = [indices; k1, k2, k3; k1, k3, k4];
            end
        end
        indices = indices + (numSlices + 1) * (numSegments + 1) * (i - 1);
        nodes.indices = [nodes.indices; indices];
    end
    
    % Save nodes to .mat file
    save('sphere_nodes.mat', 'nodes');
    
    % Plot the spheres
    % figure;
    % trisurf(nodes.indices, nodes.vertices(:, 1), nodes.vertices(:, 2), nodes.vertices(:, 3));
    % axis equal;
    % xlabel('X');
    % ylabel('Y');
    % zlabel('Z');
    % title('Spheres');
