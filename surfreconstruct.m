function depthsPX = surfreconstruct(normsPX, imgTexture)

%   input:
%       normals of all the pixels in a image, size like height x width x 3
%       and an image data for texture mapping
%   output:
%       return the depth of all pixels
%       also show the 3d surface with texture map available
%
%   description:
%       http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/#shapelet
%
%

s = size(normsPX);

slant = zeros(s(1:2));
tilt = zeros(s(1:2));

for i = 1:s(1)
    for j = 1:s(2)
        n = squeeze(normsPX(s(1)+1-i,j,:)); % take in account that image has different cooridnates system
        slant(i,j) = acos(n(3));
        tilt(i,j) = sign(n(2)) * acos(n(1));
    end
end

figure, needleplotst(slant,tilt,5,2);

depthsPX = shapeletsurf(slant, tilt, 8, 1, 2, 'slanttilt');

[x, y] = meshgrid(1:s(2), 1:s(1));
figure;

imgTexFlip = zeros(size(imgTexture));
for i = 1:3
    tmpMat = imgTexture(:,:,i);
    imgTexFlip(:,:,i) = flipud(tmpMat);     % don't know why, the texture image has to flip up-side-down
end
    
h = surf(x,y,depthsPX,'FaceColor','green','EdgeColor','none');
%set(h, 'CData', imgTexFlip, 'FaceColor', 'texturemap');
camlight left;
lighting phong;
axis vis3d;
axis off;

%figure, imshow(imgTexture), figure, imshow(imgTexFlip);
end