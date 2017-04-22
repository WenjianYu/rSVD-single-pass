%% read and process Feret photo database.
function genFeretMatrix
% function readFeret, extracting matrix from the FERET database.
% Place the program and colorferet/ in a same folder.
% Output is FeretMat.dat, a 150 GB disk file.

fp= fopen('FeretMat.dat', 'w');

processdir('colorferet/colorferet/dvd1/data/images/', fp);
disp('pass file1');
processdir('colorferet/colorferet/dvd2/data/images/', fp);

fclose(fp);
end

function processdir(dir1, matfile);
list1= ls(dir1);
list1= deblank(list1);
subdir1= regexp(list1, '\s+', 'split');
% subdir1= textscan(list1, '%s');
n1= length(subdir1);
for i=1:n1,
    dir2= subdir1{i};
    dir2= [dir1 '/' dir2];
    disp(sprintf('%s', dir2));
    list2= ls(dir2);
    %pause(0.2);
    %subdir2= textscan(list2, '%s');
    list2= deblank(list2);
    subdir2= regexp(list2, '\s+', 'split');
    n2= length(subdir2);
    for j=1:n2,
        filename= [dir2 '/' subdir2{j}];
        % process filename as a photo, and wirte to datafile
        processimg(filename, matfile);
        disp(sprintf('%d %d', i,j));
    end
end
end

function processimg(imgfile, matfile)
X= imread(imgfile);
[m,n,p]= size(X);

XX= double(X)/255;
for i=1:3,
    x= reshape(XX(:,:,i), 1, m*n);
    x= x- mean(x);
    x= x/norm(x);
    fwrite(matfile, x, 'float');
end

% random duplicate
XX= X;
for i=1:3,
    idx= randi(m*n, 1, round(m*n*0.1));
    rndvalue= randi(256, 1, length(idx))-1;
    XXi= XX(:,:, i);
    XXi(idx)= rndvalue;
    XX(:,:, i)= XXi;
end
% image(XX); axis image; axis off;   % debug
XX= double(XX)/255;
for i=1:3,
    x= reshape(XX(:,:,i), 1, m*n);
    x= x- mean(x);
    x= x/norm(x);
    fwrite(matfile, x, 'float');
end

% random duplicate
XX= X;
for i=1:3,
    idx= randi(m*n, 1, round(m*n*0.1));
    rndvalue= randi(256, 1, length(idx))-1;
    XXi= XX(:,:, i);
    XXi(idx)= rndvalue;
    XX(:,:, i)= XXi;
end
% image(XX); axis image; axis off;   % debug
XX= double(XX)/255;
for i=1:3,
    x= reshape(XX(:,:,i), 1, m*n);
    x= x- mean(x);
    x= x/norm(x);
    fwrite(matfile, x, 'float');
end
end