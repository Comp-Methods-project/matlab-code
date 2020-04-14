%%
clear
clc

a = dicomread('IM-0001-0001.dcm');
b = dicomread('IM-0001-0002.dcm');
c = dicomread('IM-0001-0003.dcm');
d = dicomread('IM-0001-0004.dcm');

figure(1)
imagesc(a)
colormap gray
axis image

figure(2)
imagesc(b)
colormap gray
axis image

figure(3)
imagesc(c)
colormap gray
axis image

figure(4)
imagesc(d)
colormap gray
axis image


%% 4/6/20 Gunjan test code.Anisotropic diffusion filter applied to figure 1
clear
clc
a = dicomread('IM-0001-0001.dcm');
b = dicomread('IM-0001-0002.dcm');
c = dicomread('IM-0001-0003.dcm');
d = dicomread('IM-0001-0004.dcm');


subplot(2,1,1)

imagesc(a)
colormap gray
axis image
% 
% figure(2)
% imagesc(b)
% colormap gray
% axis image
% 
% figure(3)
% imagesc(c)
% colormap gray
% axis image
% 
% figure(4)
% imagesc(d)
% colormap gray
% axis image

% new code from here onwards
image3d = cat(3, a,b,c,d);
% combined matrix from all 4 images. Not sure that this is all that useful
% , so far


anistropic_a=anisodiff2D(a,15,1/7,30,2); 
%function downloaded from https://www.mathworks.com/matlabcentral/fileexchange/14995-anisotropic-diffusion-perona-malik

% The idea behind image diffusion is to reduce image noise while keeping
% edges around. This will allow for easier detection of edges. 

% Might be pointless. 

subplot(2,1,2)
imagesc(anistropic_a)
colormap gray
axis image


%% 4/10/20 Gunjan 

% trying to convert anistropic_a into grayscale, but rgb2gray built in
% command is not working. Is the image already in grayscale? 



sout=imresize(anistropic_a,[256,256]);
t0=60;
th=t0+((max(anistropic_a(:))+min(anistropic_a(:)))./2);
for i=1:1:size(anistropic_a,1)
    for j=1:1:size(anistropic_a,2)
        if anistropic_a(i,j)>th
            sout(i,j)=1;
        else
            sout(i,j)=0;
        end
    end
end


%% morphological operation

% this should be extracting out the atrial fibrilation lesion area. But
% it's not. Need to thoroughly understand the code to see what's going on. 

label=bwlabel(sout);
stats=regionprops(logical(sout),'Solidity','Area','BoundingBox');
density=[stats.Solidity];
area=[stats.Area];
high_dense_area=density>0.6;
max_area=max(area(high_dense_area));
tumor_label=find(area==max_area);
tumor=ismember(label,tumor_label);

if max_area>100
   figure;
   imshow(tumor)
   title('tumor alone','FontSize',20);
else
    h = msgbox('No Tumor!!','status');
    %disp('no tumor');
    return;
end

%% 4/12/20 Andrew (updated 4/13 with imfreehand commands)

% trying out different image adjustment techniques to improve contrast

a = dicomread('IM-0001-0001.dcm');
a_imadjust = imadjust(a);
a_histeq = histeq(a);
a_adapthisteq = adapthisteq(a);

% figure(5)
% montage({a,a_imadjust,a_histeq,a_adapthisteq},'Size',[1 4])
% This was just to show the different options. Commenting out for now

figure(5)
imshow(a_imadjust) % This is the best way to improve contrast that I have found
ROI_a = imfreehand;

b_imadjust = imadjust(b);
c_imadjust = imadjust(c);
d_imadjust = imadjust(d);

figure(6)
imshow(b_imadjust)
ROI_b = imfreehand;

figure(7)
imshow(c_imadjust)
ROI_c = imfreehand;

figure(8)
imshow(d_imadjust)
ROI_d = imfreehand;

%% Tanner 4/13/20 
% this is purely for viewing the slices in 3D
D = cat(3,a,b,c,d);
figure
colormap gray
contourslice(D,[],[],[1,2,3,4,27],15);
view(3)
axis tight