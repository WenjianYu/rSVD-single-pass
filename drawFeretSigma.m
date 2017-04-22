% read singular value of Feret matrix and draw.
filename= 'FeretMat_rSVDsp_k=50_S.dat';
k= 50;
fp= fopen(filename, 'r', 'l');
sig=fread(fp, k, 'double');
fclose(fp);
plot(sig, '.-');
axis([1, 50, 0, 180]);
ylabel('\sigma_{ii}');
