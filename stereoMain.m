clear, clc;

Folder = '../data/data08';
ImgType = '*.bmp';

[Imgs, LVs] = resampling(Folder,ImgType);

%save([Folder '/resample.mat'], 'LVs', 'Imgs');

[normsPX, idxDeImg] = initialnorm(Imgs, LVs);

%save([Folder '/norms.mat'], 'normsPX');

normsOPT = refinement(normsPX);

depthsPX = surfreconstruct(normsOPT, Imgs(:,:,:,idxDeImg)/255);

%save([Folder '/depths.mat'], 'depthsPX');