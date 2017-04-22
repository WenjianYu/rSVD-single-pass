% plot eigenface from Feret.
%% read file
filename= 'FeretMat_rSVDsp_k=6_V.dat';
n= 393216;
% l= 60;
l= 6;
fp= fopen(filename, 'r', 'l');
B=fread(fp, n*l, 'double');
V= reshape(B, l, n)';
fclose(fp);


%% draw face
% i= 1;
nn= 512;
mm= 768;
subplot('position',[0,0.5,.49,.49]);
X= reshape(-V(:,1), mm, nn);
imagesc(X)
colormap gray
axis image
axis off

subplot('position',[0.5,0.5,.49,.49]);
X= reshape(V(:,2), mm, nn);
imagesc(X)
colormap gray
axis image
axis off

subplot('position',[0,0,.49,.49]);
X= reshape(V(:,3), mm, nn);
imagesc(X)
colormap gray
axis image
axis off

subplot('position',[0.5,0,.49,.49]);
X= reshape(-V(:,4), mm, nn);
imagesc(X)
colormap gray
axis image
axis off
