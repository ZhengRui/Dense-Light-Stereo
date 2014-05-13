function [RImgs, UniLV] = resampling(Folder, ImgType)

%   input:  folder name, image type, images required to be the same size
%   output: images after uniform resampling and its corresponding light
%           directions from the icosahedron vertices
%
%   description:
%           seek the nearest light di-rection L_o at one vertex in the 
%           subdivided icosahedron for each captured light direction L_i, 
%           and interpolate the image I_o at L_o by 
%
%                                                          L_o * L_i' * I_i(x,y)
%           I_o(x,y) = sum_{i | L_i's NN == L_o} * -------------------------------------
%                                                   sum_{i | L_i's NN == L_o} L_o * L_i'
%
%
%   example: resampling('../data/data02','*.bmp');

TR = IcosahedronMesh;                % generate base icosahedron
TR4 = SubdivideSphericalMesh(TR, 4); % subdivision for 4 times
IcoPts = TR4.X;                      % get coords of vertices

%sum(IcoPts(:, 3) >= 0)              % how many points has z >= 0, 
                                     % for 4th level, it should be 1313,
                                     % check with the paper
                                     
LightVec = textread([Folder '/lightvec.txt']);
LVnorm = normr(LightVec);            % normalization

NNIndex = nearestneighbour(LVnorm', IcoPts'); % find the index of nearest 
                                              % neighbour in icosahedron
                                              % vertices
[uniIndex, ~, revIndex] = unique(NNIndex);    % only part of icosahedron 
uniN = length(uniIndex);                      % vertices are the nearest
UniLV = zeros([uniN 3]);                      % neighbors, pick out them
for i = 1:length(uniIndex)
    UniLV(i,:) = IcoPts(uniIndex(i),:,:);
end
                                              
Imgs = dir([Folder '/' ImgType]);   % all images info saved in Imgs
NumImgs = size(Imgs, 1);

image = double(imread([Folder '/' Imgs(1).name]));
RImgs = zeros([size(image) uniN]);
Weight = zeros([uniN 1]);


for i = 1:NumImgs
    image = double(imread([Folder '/' Imgs(i).name]));
    %imshow(image);
    w = IcoPts(NNIndex(i),:,:) * LVnorm(i,:,:)';
    ri = revIndex(i);
    Weight(ri) = Weight(ri) + w;                    % weight based on the projection
    RImgs(:,:,:,ri) = RImgs(:,:,:,ri) + w * image;  % weight sum for each unique 
                                                    % icosahedron light direction
end

for i = 1:uniN
    RImgs(:,:,:,i) = RImgs(:,:,:,i) / Weight(i);    % weight average
end

end