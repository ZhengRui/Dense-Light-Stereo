function normsOPT = refinement(norms)

%
%   input:  initial normals after the SVD step
%   output: optimal normals (based on minimizing energy function of a 
%           markov random field, transform it into multiway cut problem)
%
%   description:  check papers referred in README of gco-3.0 library
%

s = size(norms);


TR = IcosahedronMesh;                % generate base icosahedron
TR5 = SubdivideSphericalMesh(TR, 3); % subdivision for 5 times
IcoPts = TR5.X; 

idx = IcoPts(:, 3) > 0;             % filter out light vectors with z > 0
labels = IcoPts(idx, :);
L = size(labels,1);

h = GCO_Create(s(1) * s(2), L);

dataTerm = pdist2(labels, reshape(norms, [], 3));  % reshape 2d pixel into 1d structure 
                                                   % note that matlab is column based, so the
                                                   % the 2nd pos in 1d structure is
                                                   % pixel(2,1), not pixel(1,2)

dataTerm = int32(dataTerm*10000);
%dataTerm(dataTerm > 1000) = 1000000;
                                                   
GCO_SetDataCost(h, dataTerm);

lambda = 0.5;   % weight compared to single-site term
sigma = 0.6;    % how much do you allow for the difference of neighbour normals, smaller smoother

% good parameters for each data:
%
%           lambda      sigma
%
%   data02  0.5         0.6
%   small   0.5         0.6
%   data04  0.5         0.8
%   data05  0.5         0.8
%   data06  0.5         0.7
%   data07  0.5         0.8
%   data08  0.2; 0.5    0.3; 0.4-0.6
%   data09
%   data10
%


disMat = pdist2(labels, labels);
smoothTerm = lambda * log(1 + disMat / (2 * sigma * sigma));
smoothTerm = int32(smoothTerm*10000);
%smoothTerm(smoothTerm > 1000) = 1000000;
GCO_SetSmoothCost(h, smoothTerm);

si = ones(2*s(1)*s(2)-s(1)-s(2), 1);    % following is to construct sparse neighboring matrix
sj = si;
ss = si;
k = 0;
for j = 1:s(2)
    for i = 1:s(1)
        if (i < s(1))
            k = k+1;
            si(k) = (j-1) * s(1) + i;
            sj(k) = si(k) + 1;
        end
        if (j < s(2))
            k = k+1;
            si(k) = (j-1) * s(1) + i;
            sj(k) = si(k) + s(1);
        end
    end
end

nbMat = sparse(si, sj, ss, s(1)*s(2), s(1)*s(2));
GCO_SetNeighbors(h, nbMat);

GCO_Expansion(h);
optLabels1D = GCO_GetLabeling(h);
normsOPT1D = labels(optLabels1D, :);
normsOPT = zeros(s);
for j = 1:s(2)
    for i = 1:s(1)
        normsOPT(i,j,:) = normsOPT1D((j-1)*s(1) + i, :);
    end
end


figure, imshow((-1/sqrt(3) * normsOPT(:,:,1) + 1/sqrt(3) * normsOPT(:,:,2) + 1/sqrt(3) * normsOPT(:,:,3)) / 1.1);
end