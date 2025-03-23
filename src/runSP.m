clear; clc; close all;

M = loadMesh('inputmodels\knot1.obj');
% trimesh(M.TR);
l0 = computeInitialEdgeLengths(M);

thetaT = generateTargetAngles(M);

epsilon = 1e-03;

M_ = SeamlessParametrization(M, l0, thetaT, epsilon);

visualizeParametrization(M_);



function M = loadMesh(filename)
    [vertices, faces] = readMesh(filename);
    
    M.vertices = vertices; 
    M.faces = faces;
    M.TR = triangulation(faces,vertices);
    M.edges = M.TR.edges;
    M.Fe = computeConnectivityEdges(M);
    M.edgeLengths = computeEdgeLengths(vertices, M.edges);

    M.genus = computeGenus(M);

    M.coneVertices = specifyCones(vertices, M.genus);
    M.loops = specifyDualLoops(faces, M.genus);
end

function [vertices, faces] = readMesh(filename)
    if endsWith(filename, '.obj')
        [vertices, faces] = readObj(filename);
    elseif endsWith(filename, '.stl')
        [vertices, faces] = readStl(filename);
    else
        error('Unsupported file format. Use .obj or .stl.');
    end
end

function [vertices, faces] = readObj(filename)
    fid = fopen(filename, 'r');
    vertices = []; faces = [];

    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(line)
            parts = strsplit(line, ' ');
            if parts{1} == 'v' 
                vertices = [vertices; str2double(parts(2:end))];
            elseif parts{1} == 'f'  
                faces = [faces; str2double(parts(2:end))]; 
            end
        end
    end
    fclose(fid);
end

function [vertices, faces] = readStl(filename)
    fid = fopen(filename, 'r');
    vertices = []; faces = [];
    
    fclose(fid);
end

function Fe = computeConnectivityEdges(M)
    E = M.edges; NE = size(E,1);
    F = M.faces; 
    
    Fe = int32(zeros(size(F)));
    for ei = 1:NE
        fidxes = edgeAttachments(M.TR,E(ei,1),E(ei,2));
        fi1 = fidxes{:}(1); fi2 = fidxes{:}(2);

        Fe(fi1,F(fi1,:) ~= E(ei,1) & F(fi1,:) ~= E(ei,2)) = ei;
        Fe(fi2,F(fi2,:) ~= E(ei,1) & F(fi2,:) ~= E(ei,2)) = ei;
    end
end

function edgeLengths = computeEdgeLengths(vertices, edges)
    NE = size(edges, 1);
    edgeLengths = zeros(NE, 1);
    for i = 1:NE
        v1 = vertices(edges(i, 1), :);
        v2 = vertices(edges(i, 2), :);
        edgeLengths(i) = norm(v1 - v2);
    end
end

function genus = computeGenus(M)
    % 使用欧拉公式：V - E + F = 2 - 2g
    V = size(M.vertices, 1); 
    E = size(M.edges, 1);   
    F = length(M.faces);    
    
    genus = (2 - (V - E + F)) / 2;
end

function coneVertices = specifyCones(vertices, genus)
    NV = size(vertices, 1);
    
    if genus == 0
        numCones = 2;
    else
        numCones = 2 * genus - 2;
    end

    coneVertices = randperm(NV, numCones); 
end

function loops = specifyDualLoops(faces, genus)
    NL = 2 * genus;
    loops = cell(1, NL);

    dualGraph = buildDualGraph(faces);
    hasvisited = false(size(faces, 1), 1);
   
    for i = 1:NL
        unvisitedFaces = find(~hasvisited);
        while true 
            startFace = unvisitedFaces(randi(length(unvisitedFaces)));
            neighbors = dualGraph(startFace, :);
            neighbors = neighbors(~hasvisited(neighbors));
            if length(neighbors) == 3
                break;
            end
        end
         
        loop = dfsDualLoop(dualGraph, startFace, hasvisited);
        loops{i} = loop;
    end
end

function dualGraph = buildDualGraph(faces)
    NF = size(faces, 1);
    dualGraph = zeros(NF, 3);

    for fi = 1:NF
        neighbors = [];
        for fj = 1:NF
            if fi ~= fj && shareEdge(faces(fi,:), faces(fj,:))
                neighbors = [neighbors, fj];
            end
        end
        dualGraph(fi,:) = neighbors;
    end
end

function shared = shareEdge(face1, face2)
    shared = false;
    for i = 1:3
        for j = 1:3
            if face1(i) == face2(j) && face1(mod(i,3)+1) == face2(mod(j-2,3)+1)
                shared = true;
                return;
            end
        end
    end
end

function [loop, hasvisited] = dfsDualLoop(dualGraph, startFace, hasvisited)
    visited = false(size(dualGraph, 1), 1);
    stack = startFace; 
    loop = []; 
    
    while ~isempty(stack)
        currentFace = stack(end);
        stack(end) = [];
        
        visited(currentFace) = true; 
        loop = [loop, currentFace]; 
        
        neighbors = dualGraph(currentFace, :);
        neighbors = neighbors(~visited(neighbors) & ~hasvisited(neighbors));
        
        if length(neighbors) >= 2
            stack = [stack, neighbors(randi(length(neighbors)))];
        else
            if ~isempty(loop) && ismember(startFace,dualGraph(loop(end), :))
                hasvisited(loop) = true;
            else
                loop = []; 
                visited(:) = false; 
                stack = startFace;
            end
        end
    end
end

function l0 = computeInitialEdgeLengths(M)
    l0 = M.edgeLengths;
end

function thetaT = generateTargetAngles(M)
    NV = size(M.vertices, 1);
    thetaT.vertexAngles = 2 * pi * ones(NV - 1, 1);
    coneIndices = M.coneVertices;

    for i = 1:length(coneIndices)
        if mod(i, 2) == 1
            thetaT.vertexAngles(coneIndices(i)) = 3*pi/2;
        else
            thetaT.vertexAngles(coneIndices(i)) = 5*pi/2;
        end
    end
    
    thetaT.loopAngles = (2 * pi) * ones(2 * M.genus, 1);
end

function visualizeParametrization(M_)
    
end