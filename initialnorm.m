function [normsPX, idxDeImg] = initialnorm(Imgs, LVs)

%   input:  a big images matrix after uniform resampling, dimension should
%           be [size(img) #images], Imgs(:,:,:,i) saves ith image; LVs are the
%           corresponding icosahedra vertices, which are light directions
%   output: initial norm direction after SVD, and the index of denominator
%           image
%
%   description:
%           Lambertian model: 
%
%                    I_i     N * L_i
%                   ----- = ---------
%                    I_d     N * L_d
%
%           can be transform into form of A_x*N_x + A_y*N_y + A_z*N_z = 0
%           take one image as reference, which is called denominator image,
%           construct k-1 equation, and solve overdeterminant equations
%
%                        AX = 0
%
%            using SVD, A = U*S*V', S save singular values, while the
%            solution of X is the last column of V
%
%           Criteria for choosing denominator image: for each resampled image_i, count the number of pixels whose intensity rank
%           satisfies rank_i > L, where L ≥ median. Let k_i be the total number of pixels in image_i 
%           satisfying this condition, r_i be the mean rank among the pixels in image_i that satisfies
%           this condition. The denominator image is defined to be one
%           with 1) maximum k_i and 2) r_i lower than some thresh H. 
%           Currently, I set L and H to be the 70th and 90th percentiles
%           respectively.
%   
%

s = size(Imgs);
grayImgs = zeros([s(1:2) s(4)]);

% to get intensity of each pixel, first convert to gray images
for i = 1:s(4)
    grayImgs(:,:,i) = 0.2989 * Imgs(:,:,1,i) + 0.5870 * Imgs(:,:,2,i) + 0.1140 * Imgs(:,:,3,i);
%    figure, imshow(Imgs(:,:,:,i)/255), figure, imshow(grayImgs(:,:,:,i)/255); % check if convertion correct
end

rankPixel = zeros([s(1:2) s(4)]);
rankIJ = zeros([s(4) 1]);

for i = 1:s(1)
    for j = 1:s(2)
        pxIJ = grayImgs(i,j,:);             % take intensity of pixel(i,j) from all images
        [pxIJSort, ind] = sort(pxIJ);       % sort it, get the sort index, pxIJSort = pxIJ(ind)
        rankIJ(ind) = 1:s(4);               % reverse sort, pxIJ = pxIJSort(rankIJ)
        rankPixel(i,j,:) = rankIJ;
    end
end

% select the denominator image
L = 0.7 * s(4);     % pixel intensity rank lower bound
H = 0.9 * s(4);     % average intensity rank upper bound
stat = zeros([s(4) 2]);

for i = 1:s(4)
    effInd = rankPixel(:,:,i) >= L;         % take out index of pixels with rank > L
    stat(i,1) = sum(effInd(:));             % count number
    effPixel = rankPixel(:,:,i) .* effInd;  % take out pixels with rank > L
    stat(i,2) = sum(effPixel(:)) / stat(i,1);           % get average intensity rank
end

[val idxDeImg] = max(stat(:,1) .* (stat(:,2) < H)); % find the denominator image index

figure, imshow(Imgs(:,:,:,idxDeImg)/255);
figure, imshow(grayImgs(:,:,idxDeImg)/255); % view the denominator image


lv_de = LVs(idxDeImg, :);
lv_ks = [LVs(1:(idxDeImg-1),:); LVs((idxDeImg+1):end,:)];
normsPX = zeros([s(1:2) 3]);
for i = 1:s(1)
    for j = 1:s(2)
        pxIJ = squeeze(grayImgs(i,j,:));
        I_de = pxIJ(idxDeImg);
        I_ks = [pxIJ(1:(idxDeImg-1)); pxIJ((idxDeImg+1):end)];
        A = I_de * lv_ks - I_ks * lv_de;  % A(k) = I(i,j,d) * L(k) - I(i,j,k)*L(d)
        [~,~,v] = svd(A,0);               % svd, ||N|| == 1 satisfied naturally
        sgn = sign(v(3,3));               % z component maybe negative
        normsPX(i,j,:) = sgn * v(:,3);
    end
end

figure, imshow((-1/sqrt(3) * normsPX(:,:,1) + 1/sqrt(3) * normsPX(:,:,2) + 1/sqrt(3) * normsPX(:,:,3)) / 1.1);


end