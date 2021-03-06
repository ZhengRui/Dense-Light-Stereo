clear,clc;

p1 = [zeros(1,15) [0:.2:10] zeros(1,15)];
p2 = zeros(1,length(p1));
n = ceil(length(p1)/2);
p3 = [1:n n:-1:1]/n*6;
p3 = p3(1:length(p1));
    
z = [ones(15,1)*p2
     ones(15,1)*p1
     ones(15,1)*p3
     ones(15,1)*p2];
 
[dx, dy] = gradient(z);
[s,t] = grad2slanttilt(dx,dy);

figure(1); clf;% surfl(z), shading interp; colormap(copper);
surf(z);       % colormap(white);
figure(2),needleplotst(s,t,5,2), axis('off');

depthsPX = shapeletsurf(s, t, 7, 1, 2, 'slanttilt');
sz=size(s);
[x, y] = meshgrid(1:sz(2), 1:sz(1));
figure(3);
surf(x,y,depthsPX,'FaceColor','green','EdgeColor','none');
camlight left; lighting phong;