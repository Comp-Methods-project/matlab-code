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


%% Gunjan 4/14/20 
% Code that outputs mean value of ROI, SD of ROI, area in pixels,
% perimeter, Centroid, and Center of mass


clear
clc
% building on top of Andrew's brightness image code

a = dicomread('IM-0001-0001.dcm');
a_imadjust = imadjust(a);
a_histeq = histeq(a);
a_adapthisteq = adapthisteq(a);



figure(5)
imshow(a_imadjust)% This is the best way to improve contrast that I have found
% NEW CODE FROM HERE ONWARDS
axis on;
title('Original MRI image');
set(gcf,'Position',get(0,'Screensize')); %Maximizing the figure

%Now asking user to draw on the image, after which we will apply a mask
message =sprintf('Left click the mouse to begin drawing.\n Stop holding the mouse button to finish');
uiwait(msgbox(message));
ROI_a = imfreehand; %Line of code from Andrew that does the drawing
binaryImage = ROI_a.createMask; %This is creating a mask from the ROI
%Mask is basically a binary image of the ROI
xy=ROI_a.getPosition;


%using subplots to show more images
subplot(2,3,1);
imshow(a_imadjust,[]);
axis on;
drawnow;
title('Original MRI image');


% Now we display the mask that drawn by user
subplot(2,3,2);
imshow(binaryImage);
axis on;
title('Binary mask of the MRI image');

% Time to label the MRI image and compute the centroid and center of mass
labeledImage=bwlabel(binaryImage);
measurements = regionprops(binaryImage,a_imadjust, 'area','Centroid','WeightedCentroid','Perimeter');
area = measurements.Area;
centroid = measurements.Centroid;
centerofMass = measurements.WeightedCentroid;
perimeter = measurements.Perimeter;


% Calculating the area from the pixels selected by user
numberofPixels1=sum(binaryImage(:))
numberofPixels2=bwarea(binaryImage)

% Getting the coordinates of the boundary of MRI region selected by the
% user
structBoundaries = bwboundaries(binaryImage);
xy= structBoundaries{1}; % This gives us a n by 2 arry of x and y coordinates
x = xy(:,2); % Column of interest
y = xy(:,1); % row of interest
subplot(2,3,1); % plotting over the original image
hold on;
plot(x,y,'LineWidth',2);
drawnow; % forcing matlab to draw it asap


% Now we are burning line into image by setting it to 255 wherever the mask
% is true
burnedImage = a_imadjust;
burnedImage(binaryImage) = 255;


% Here we are displaying the image with the mask burned in
subplot(2,3,3);
imshow(burnedImage);
axis on;
caption = sprintf('New image with mask \n burned into image');
title(caption);


%Now keeping only the part of the image that is inside the mask
blackMaskedImage=a_imadjust;
blackMaskedImage(~binaryImage) = 0;
subplot(2,3,4);
imshow(blackMaskedImage);
axis on;
title('Masked Outside Region');

% Calculating the mean
meanGL=mean(blackMaskedImage(binaryImage));
sdGL = std(double(blackMaskedImage(binaryImage)));

% Placing markeers at the centroid and center of mass
hold on;
plot(centroid(1),centroid(2),'b+','MarkerSize',20,'LineWidth',2);
plot(centerofMass(1),centerofMass(2),'g+','MarkerSize',10,'LineWidth',2);

% Blackening inside the region
insideMasked = a_imadjust;
insideMasked(binaryImage) = 0;
subplot(2,3,5);
imshow(insideMasked);
axis on;
title('Masked inside region');


% Cropping the image
leftColumn = min(x);
rightColumn = max(x);
topLine = min(y);
bottomLine = max(y);
width = rightColumn - leftColumn + 1;
height = bottomLine - topLine + 1;
croppedImage = imcrop(blackMaskedImage,[leftColumn,topLine,width,height]);


% displaying cropped image
subplot(2,3,6);
imshow(croppedImage);
axis 'on'
title('Cropped Image');

% Placing crosses at the centroid and center of mass
hold on;
plot(centroid(1)-leftColumn,centroid(2)-topLine,'b+','MarkerSize',20,'LineWidth',2);
plot(centerofMass(1)-leftColumn,centerofMass(2)-topLine,'g+','MarkerSize',20,'LineWidth',2);

% report the results of the calculation
message=sprintf('Mean value of ROI = %.3f\n SD of ROI = %.3f\nNumber of pixels =%d\nArea in pixels=%.2f\nPerimeter = %.2f\nCentroid @ (x,y) = (%.1f,%.1f)\n Center of Mass @ (x,y) = (%.1f,%.1f)\nBlue crosshairs @ centroid.\n Green crosshairs @ center of mass.', meanGL,sdGL,numberofPixels1,numberofPixels2,perimeter,centroid(1),centroid(2),centerofMass(1),centerofMass(2));
msgbox(message);





% b_imadjust = imadjust(b);
% c_imadjust = imadjust(c);
% d_imadjust = imadjust(d);
% 
% figure(6)
% imshow(b_imadjust)
% ROI_b = drawfreehand;
% 
% figure(7)
% imshow(c_imadjust)
% ROI_c = drawfreehand;
% 
% figure(8)
% imshow(d_imadjust)
% ROI_d = drawfreehand;

