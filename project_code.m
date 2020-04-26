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

info1 = dicominfo('IM-0001-0001.dcm'); % this is to extract metadata about the width and height of the images
info2 = dicominfo('IM-0001-0002.dcm');
info3 = dicominfo('IM-0001-0003.dcm');
info4=  dicominfo('IM-0001-0004.dcm');



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
imshow(a_imadjust)
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


% Calculating the area and volume from the pixels selected by user
numberofPixels1=sum(binaryImage(:))
numberofPixels2=bwarea(binaryImage)
Volume = (.625^2)*4*numberofPixels1

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
message=sprintf('Mean value of ROI = %.3f\n SD of ROI = %.3f\nNumber of pixels =%d\nVolume of ROI=%.2f', meanGL,sdGL,numberofPixels1,Volume);
msgbox(message);

%% so here is the same code as above repeated for each image. I know you guys are using DICOM so just switch that out for the png.
b=imread('2.png');
b_imadjust = imadjust(b);
b_histeq = histeq(b);
b_adapthisteq = adapthisteq(b);
figure(2)
imshow(b_imadjust)% This is the best way to improve contrast that I have found
% NEW CODE FROM HERE ONWARDS
axis on;
title('Image B');
set(gcf,'Position',get(0,'Screensize')); %Maximizing the figure
%Now asking user to draw on the image, after which we will apply a mask
message =sprintf('Left click the mouse to begin drawing.\n Stop holding the mouse button to finish');
uiwait(msgbox(message));
ROI_b = imfreehand; %Line of code from Andrew that does the drawing
binaryImage = ROI_b.createMask; %This is creating a mask from the ROI
%Mask is basically a binary image of the ROI
xy=ROI_b.getPosition;
%using subplots to show more images
subplot(2,3,1);
imshow(b_imadjust,[]);
axis on;
drawnow;
title('Image B');
% Now we display the mask that drawn by user
subplot(2,3,2);
imshow(binaryImage);
axis on;
title('Binary mask of the MRI image');
% Time to label the MRI image and compute the centroid and center of mass
labeledImage=bwlabel(binaryImage);
measurements = regionprops(binaryImage,b_imadjust, 'area','Centroid','WeightedCentroid','Perimeter');
area = measurements.Area;
centroid = measurements.Centroid;
centerofMass = measurements.WeightedCentroid;
perimeter = measurements.Perimeter;
% Calculating the area from the pixels selected by user
numberofPixels1=sum(binaryImage(:))
numberofPixels2b=bwarea(binaryImage)
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
burnedImage = b_imadjust;
burnedImage(binaryImage) = 255;
% Here we are displaying the image with the mask burned in
subplot(2,3,3);
imshow(burnedImage);
axis on;
caption = sprintf('New image with mask \n burned into image');
title(caption);
%Now keeping only the part of the image that is inside the mask
blackMaskedImage=b_imadjust;
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
insideMasked = b_imadjust;
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
message=sprintf('Mean value of ROI = %.3f\n SD of ROI = %.3f\nNumber of pixels =%d\nArea in pixels=%.2f\nPerimeter = %.2f\nCentroid @ (x,y) = (%.1f,%.1f)\n Center of Mass @ (x,y) = (%.1f,%.1f)\nBlue crosshairs @ centroid.\n Green crosshairs @ center of mass.', meanGL,sdGL,numberofPixels1,numberofPixels2b,perimeter,centroid(1),centroid(2),centerofMass(1),centerofMass(2));
msgbox(message);


c=imread('3.png');
c_imadjust = imadjust(c);
c_histeq = histeq(c);
c_adapthisteq = adapthisteq(c);
figure(3)
imshow(c_imadjust)
axis on;
title('Image C');
set(gcf,'Position',get(0,'Screensize')); %Maximizing the figure
message =sprintf('Left click the mouse to begin drawing.\n Stop holding the mouse button to finish');
uiwait(msgbox(message));
ROI_c = imfreehand; %Line of code from Andrew that does the drawing
binaryImage = ROI_c.createMask; %This is creating a mask from the ROI
xy=ROI_c.getPosition;
subplot(2,3,1);
imshow(c_imadjust,[]);
axis on;
drawnow;
title('Image C');
subplot(2,3,2);
imshow(binaryImage);
axis on;
title('Binary mask of the MRI image');
labeledImage=bwlabel(binaryImage);
measurements = regionprops(binaryImage,c_imadjust, 'area','Centroid','WeightedCentroid','Perimeter');
area = measurements.Area;
centroid = measurements.Centroid;
centerofMass = measurements.WeightedCentroid;
perimeter = measurements.Perimeter;
numberofPixels1=sum(binaryImage(:))
numberofPixels2c=bwarea(binaryImage)
structBoundaries = bwboundaries(binaryImage);
xy= structBoundaries{1}; % This gives us a n by 2 arry of x and y coordinates
x = xy(:,2); % Column of interest
y = xy(:,1); % row of interest
subplot(2,3,1); % plotting over the original image
hold on;
plot(x,y,'LineWidth',2);
drawnow; % forcing matlab to draw it asap
burnedImage = c_imadjust;
burnedImage(binaryImage) = 255;
subplot(2,3,3);
imshow(burnedImage);
axis on;
caption = sprintf('New image with mask \n burned into image');
title(caption);
blackMaskedImage=c_imadjust;
blackMaskedImage(~binaryImage) = 0;
subplot(2,3,4);
imshow(blackMaskedImage);
axis on;
title('Masked Outside Region');
meanGL=mean(blackMaskedImage(binaryImage));
sdGL = std(double(blackMaskedImage(binaryImage)));
hold on;
plot(centroid(1),centroid(2),'b+','MarkerSize',20,'LineWidth',2);
plot(centerofMass(1),centerofMass(2),'g+','MarkerSize',10,'LineWidth',2);
insideMasked = c_imadjust;
insideMasked(binaryImage) = 0;
subplot(2,3,5);
imshow(insideMasked);
axis on;
title('Masked inside region');
leftColumn = min(x);
rightColumn = max(x);
topLine = min(y);
bottomLine = max(y);
width = rightColumn - leftColumn + 1;
height = bottomLine - topLine + 1;
croppedImage = imcrop(blackMaskedImage,[leftColumn,topLine,width,height]);
subplot(2,3,6);
imshow(croppedImage);
axis 'on'
title('Cropped Image');
hold on;
plot(centroid(1)-leftColumn,centroid(2)-topLine,'b+','MarkerSize',20,'LineWidth',2);
plot(centerofMass(1)-leftColumn,centerofMass(2)-topLine,'g+','MarkerSize',20,'LineWidth',2);
message=sprintf('Mean value of ROI = %.3f\n SD of ROI = %.3f\nNumber of pixels =%d\nArea in pixels=%.2f\nPerimeter = %.2f\nCentroid @ (x,y) = (%.1f,%.1f)\n Center of Mass @ (x,y) = (%.1f,%.1f)\nBlue crosshairs @ centroid.\n Green crosshairs @ center of mass.', meanGL,sdGL,numberofPixels1,numberofPixels2c,perimeter,centroid(1),centroid(2),centerofMass(1),centerofMass(2));
msgbox(message);


d=imread('4.png');
d_imadjust = imadjust(d);
d_histeq = histeq(d);
a_adapthisteq = adapthisteq(d);
figure(4)
imshow(d_imadjust)
axis on;
title('Image D');
set(gcf,'Position',get(0,'Screensize')); %Maximizing the figure
message =sprintf('Left click the mouse to begin drawing.\n Stop holding the mouse button to finish');
uiwait(msgbox(message));
ROI_d = imfreehand; %Line of code from Andrew that does the drawing
binaryImage = ROI_d.createMask; %This is creating a mask from the ROI
xy=ROI_d.getPosition;
subplot(2,3,1);
imshow(d_imadjust,[]);
axis on;
drawnow;
title('Image D');
subplot(2,3,2);
imshow(binaryImage);
axis on;
title('Binary mask of the MRI image');
labeledImage=bwlabel(binaryImage);
measurements = regionprops(binaryImage,d_imadjust, 'area','Centroid','WeightedCentroid','Perimeter');
area = measurements.Area;
centroid = measurements.Centroid;
centerofMass = measurements.WeightedCentroid;
perimeter = measurements.Perimeter;
numberofPixels1=sum(binaryImage(:))
numberofPixels2d=bwarea(binaryImage)
structBoundaries = bwboundaries(binaryImage);
xy= structBoundaries{1}; % This gives us a n by 2 arry of x and y coordinates
x = xy(:,2); % Column of interest
y = xy(:,1); % row of interest
subplot(2,3,1); % plotting over the original image
hold on;
plot(x,y,'LineWidth',2);
drawnow; % forcing matlab to draw it asap
burnedImage = d_imadjust;
burnedImage(binaryImage) = 255;
subplot(2,3,3);
imshow(burnedImage);
axis on;
caption = sprintf('New image with mask \n burned into image');
title(caption);
blackMaskedImage=d_imadjust;
blackMaskedImage(~binaryImage) = 0;
subplot(2,3,4);
imshow(blackMaskedImage);
axis on;
title('Masked Outside Region');
meanGL=mean(blackMaskedImage(binaryImage));
sdGL = std(double(blackMaskedImage(binaryImage)));
hold on;
plot(centroid(1),centroid(2),'b+','MarkerSize',20,'LineWidth',2);
plot(centerofMass(1),centerofMass(2),'g+','MarkerSize',10,'LineWidth',2);
insideMasked = d_imadjust;
insideMasked(binaryImage) = 0;
subplot(2,3,5);
imshow(insideMasked);
axis on;
title('Masked inside region');
leftColumn = min(x);
rightColumn = max(x);
topLine = min(y);
bottomLine = max(y);
width = rightColumn - leftColumn + 1;
height = bottomLine - topLine + 1;
croppedImage = imcrop(blackMaskedImage,[leftColumn,topLine,width,height]);
subplot(2,3,6);
imshow(croppedImage);
axis 'on'
title('Cropped Image');
hold on;
plot(centroid(1)-leftColumn,centroid(2)-topLine,'b+','MarkerSize',20,'LineWidth',2);
plot(centerofMass(1)-leftColumn,centerofMass(2)-topLine,'g+','MarkerSize',20,'LineWidth',2);
message=sprintf('Mean value of ROI = %.3f\n SD of ROI = %.3f\nNumber of pixels =%d\nArea in pixels=%.2f\nPerimeter = %.2f\nCentroid @ (x,y) = (%.1f,%.1f)\n Center of Mass @ (x,y) = (%.1f,%.1f)\nBlue crosshairs @ centroid.\n Green crosshairs @ center of mass.', meanGL,sdGL,numberofPixels1,numberofPixels2d,perimeter,centroid(1),centroid(2),centerofMass(1),centerofMass(2));
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



%% 4/24/20 Gunjan update

% run the top segment first; and then initialize below


D = cat(3,a,b,c,d);
cubic_pixels = nnz(D)
width = 512; % width and height of all 4 images are the same
height = 368;
bitDepth = 12;
% Bit detph represents number of integer levels used to encode intensity in
% an image. A 1 bit image has only 2 possible values: 0 or 1. 2 bit image
% has 4 possible values: 0,1,2, and 3. With 12 bit image, we have 4096
% possible values. 

sliceThickness = 4; %in mm
sliceSpacing = 4; % in mm
% slice spacing is called "spacing between splices" in dicom metadata.
% After lots of googling, it is my understanding that when spacing of
% splices and slicethickness are the same value, it means that the slices
% were taken contiguously. As in, there is no gap between the slices. If
% the slice spacing value was smaller than thickness, that would have meant
% overlap between slices. If the spacing value was larger, that would have
% meant gap between slices. 

% We only need to consider splice thickness for z value

info2.PixelSpacing;
% ans =
% 
%     0.6250 ; in mm
%     0.6250 ; in mm


% PixelSpacing is another parameter available in dicom metadata.
% The first value is the vertical spacing between two adjacent pixels
% The second value is horizontal spacing between two adjacent pixels.

% pixelspacing squared * 4mm  = mm cube /pixel 

%
%g input in a while loop



%% this is what I used as the callback for the image slider. just add a new slider and go to the code view and you will see something like 
% value = app.ImageSelectorSlider.Value; and just copy and paste the code below and the slider SHOULD work. The code for the other buttons
% you should already have. 
a = dicomread('IM-0001-0001.dcm');
b = dicomread('IM-0001-0002.dcm');
c = dicomread('IM-0001-0003.dcm');
d = dicomread('IM-0001-0004.dcm');
if value>.5 & value<1.5;
    a_imadjust = imadjust(a);
    a_histeq = histeq(a);
    a_adapthisteq = adapthisteq(a);    
    imshow(a_imadjust);
    
    axis on;
    title('Image 1');
    set(gcf,'Position',get(0,'Screensize')); %Maximizing the figure
    message =sprintf('(Image 1) Left click the mouse to begin drawing.\n Stop holding the mouse button to finish');
    uiwait(msgbox(message));
    ROI_a = imfreehand; %Line of code from Andrew that does the drawing
    binaryImage = ROI_a.createMask; %This is creating a mask from the ROI
    xy=ROI_a.getPosition;
    subplot(2,3,1);
    imshow(a_imadjust,[]);
    axis on;
    drawnow;
    title('Image 1');
    subplot(2,3,2);
    imshow(binaryImage);
    axis on;
    title('Binary mask of the Heart image');
    labeledImage=bwlabel(binaryImage);
    measurements = regionprops(binaryImage,a_imadjust, 'area','Centroid','WeightedCentroid','Perimeter');
    area = measurements.Area;
    centroid = measurements.Centroid;
    centerofMass = measurements.WeightedCentroid;
    perimeter = measurements.Perimeter;
    numberofPixels1=sum(binaryImage(:))
    numberofPixels2=bwarea(binaryImage)
    structBoundaries = bwboundaries(binaryImage);
    xy= structBoundaries{1}; % This gives us a n by 2 arry of x and y coordinates
    x = xy(:,2); % Column of interest
    y = xy(:,1); % row of interest
    subplot(2,3,1); % plotting over the original image
    hold on;
    plot(x,y,'LineWidth',2);
    drawnow; % forcing matlab to draw it asap
    burnedImage = a_imadjust;
    burnedImage(binaryImage) = 255;
    subplot(2,3,3);
    imshow(burnedImage);
    axis on;
    caption = sprintf('New image with mask \n burned into image');
    title(caption);
    blackMaskedImage=a_imadjust;
    blackMaskedImage(~binaryImage) = 0;
    subplot(2,3,4);
    imshow(blackMaskedImage);
    axis on;
    title('Masked Outside Region');
    meanGL=mean(blackMaskedImage(binaryImage));
    sdGL = std(double(blackMaskedImage(binaryImage)));
    hold on;
    plot(centroid(1),centroid(2),'b+','MarkerSize',20,'LineWidth',2);
    plot(centerofMass(1),centerofMass(2),'g+','MarkerSize',10,'LineWidth',2);
    insideMasked = a_imadjust;
    insideMasked(binaryImage) = 0;
    subplot(2,3,5);
    imshow(insideMasked);
    axis on;
    title('Masked inside region');
    leftColumn = min(x);
    rightColumn = max(x);
    topLine = min(y);
    bottomLine = max(y);
    width = rightColumn - leftColumn + 1;
    height = bottomLine - topLine + 1;
    croppedImage = imcrop(blackMaskedImage,[leftColumn,topLine,width,height]);
    subplot(2,3,6);
    imshow(croppedImage);
    axis 'on'
    title('Cropped Image');
    hold on;
    plot(centroid(1)-leftColumn,centroid(2)-topLine,'b+','MarkerSize',20,'LineWidth',2);
    plot(centerofMass(1)-leftColumn,centerofMass(2)-topLine,'g+','MarkerSize',20,'LineWidth',2);
    volume=(.625^2)*4*numberofPixels1;
    message=sprintf('Mean value of ROI = %.3f\n SD of ROI = %.3f\nNumber of pixels =%d\nArea in pixels=%.2f\nPerimeter = %.2f\nCentroid @ (x,y) = (%.1f,%.1f)\n Center of Mass @ (x,y) = (%.1f,%.1f)\nBlue crosshairs @ centroid.\n Green crosshairs @ center of mass.\nVolume=%.2fmm^3', meanGL,sdGL,numberofPixels1,numberofPixels2,perimeter,centroid(1),centroid(2),centerofMass(1),centerofMass(2),volume);
    msgbox(message);

elseif value>1.5 & value<2.5;
    b_imadjust = imadjust(b);
    b_histeq = histeq(b);
    b_adapthisteq = adapthisteq(b);    
    imshow(b_imadjust);
    
    axis on;
    title('Image 2');
    set(gcf,'Position',get(0,'Screensize')); %Maximizing the figure
    %Now asking user to draw on the image, after which we will apply a mask
    message =sprintf('(Image 2) Left click the mouse to begin drawing.\n Stop holding the mouse button to finish');
    uiwait(msgbox(message));
    ROI_b = imfreehand; %Line of code from Andrew that does the drawing
    binaryImage = ROI_b.createMask; %This is creating a mask from the ROI
    %Mask is basically a binary image of the ROI
    xy=ROI_b.getPosition;
    %using subplots to show more images
    subplot(2,3,1);
    imshow(b_imadjust,[]);
    axis on;
    drawnow;
    title('Image 2');
    % Now we display the mask that drawn by user
    subplot(2,3,2);
    imshow(binaryImage);
    axis on;
    title('Binary mask of the Heart image');
    % Time to label the MRI image and compute the centroid and center of mass
    labeledImage=bwlabel(binaryImage);
    measurements = regionprops(binaryImage,b_imadjust, 'area','Centroid','WeightedCentroid','Perimeter');
    area = measurements.Area;
    centroid = measurements.Centroid;
    centerofMass = measurements.WeightedCentroid;
    perimeter = measurements.Perimeter;
    % Calculating the area from the pixels selected by user
    numberofPixels1=sum(binaryImage(:))
    numberofPixels2b=bwarea(binaryImage)
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
    burnedImage = b_imadjust;
    burnedImage(binaryImage) = 255;
    % Here we are displaying the image with the mask burned in
    subplot(2,3,3);
    imshow(burnedImage);
    axis on;
    caption = sprintf('New image with mask \n burned into image');
    title(caption);
    %Now keeping only the part of the image that is inside the mask
    blackMaskedImage=b_imadjust;
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
    insideMasked = b_imadjust;
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
    volume=(.625^2)*4*numberofPixels1;
    message=sprintf('Mean value of ROI = %.3f\n SD of ROI = %.3f\nNumber of pixels =%d\nArea in pixels=%.2f\nPerimeter = %.2f\nCentroid @ (x,y) = (%.1f,%.1f)\n Center of Mass @ (x,y) = (%.1f,%.1f)\nBlue crosshairs @ centroid.\n Green crosshairs @ center of mass.\nVolume=%.2fmm^3', meanGL,sdGL,numberofPixels1,numberofPixels2b,perimeter,centroid(1),centroid(2),centerofMass(1),centerofMass(2),volume);
    msgbox(message);
    
elseif value>2.5 & value<3.5;
    c_imadjust = imadjust(c);
    c_histeq = histeq(c);
    c_adapthisteq = adapthisteq(c);
    imshow(c_imadjust);
    title('Image 3');
    set(gcf,'Position',get(0,'Screensize')); %Maximizing the figure
    message =sprintf('(Image 3) Left click the mouse to begin drawing.\n Stop holding the mouse button to finish');
    uiwait(msgbox(message));
    ROI_c = imfreehand; %Line of code from Andrew that does the drawing
    binaryImage = ROI_c.createMask; %This is creating a mask from the ROI
    xy=ROI_c.getPosition;
    subplot(2,3,1);
    imshow(c_imadjust,[]);
    axis on;
    drawnow;
    title('Image 3');
    subplot(2,3,2);
    imshow(binaryImage);
    
    axis on;
    title('Binary mask of the Heart image');
    labeledImage=bwlabel(binaryImage);
    measurements = regionprops(binaryImage,c_imadjust, 'area','Centroid','WeightedCentroid','Perimeter');
    area = measurements.Area;
    centroid = measurements.Centroid;
    centerofMass = measurements.WeightedCentroid;
    perimeter = measurements.Perimeter;
    numberofPixels1=sum(binaryImage(:))
    numberofPixels2c=bwarea(binaryImage)
    structBoundaries = bwboundaries(binaryImage);
    xy= structBoundaries{1}; % This gives us a n by 2 arry of x and y coordinates
    x = xy(:,2); % Column of interest
    y = xy(:,1); % row of interest
    subplot(2,3,1); % plotting over the original image
    hold on;
    plot(x,y,'LineWidth',2);
    drawnow; % forcing matlab to draw it asap
    burnedImage = c_imadjust;
    burnedImage(binaryImage) = 255;
    subplot(2,3,3);
    imshow(burnedImage);
    axis on;
    caption = sprintf('New image with mask \n burned into image');
    title(caption);
    blackMaskedImage=c_imadjust;
    blackMaskedImage(~binaryImage) = 0;
    subplot(2,3,4);
    imshow(blackMaskedImage);
    axis on;
    title('Masked Outside Region');
    meanGL=mean(blackMaskedImage(binaryImage));
    sdGL = std(double(blackMaskedImage(binaryImage)));
    hold on;
    plot(centroid(1),centroid(2),'b+','MarkerSize',20,'LineWidth',2);
    plot(centerofMass(1),centerofMass(2),'g+','MarkerSize',10,'LineWidth',2);
    insideMasked = c_imadjust;
    insideMasked(binaryImage) = 0;
    subplot(2,3,5);
    imshow(insideMasked);
    axis on;
    title('Masked inside region');
    leftColumn = min(x);
    rightColumn = max(x);
    topLine = min(y);
    bottomLine = max(y);
    width = rightColumn - leftColumn + 1;
    height = bottomLine - topLine + 1;
    croppedImage = imcrop(blackMaskedImage,[leftColumn,topLine,width,height]);
    subplot(2,3,6);
    imshow(croppedImage);
    axis 'on'
    title('Cropped Image');
    hold on;
    plot(centroid(1)-leftColumn,centroid(2)-topLine,'b+','MarkerSize',20,'LineWidth',2);
    plot(centerofMass(1)-leftColumn,centerofMass(2)-topLine,'g+','MarkerSize',20,'LineWidth',2);
    volume=(.625^2)*4*numberofPixels1;
    message=sprintf('Mean value of ROI = %.3f\n SD of ROI = %.3f\nNumber of pixels =%d\nArea in pixels=%.2f\nPerimeter = %.2f\nCentroid @ (x,y) = (%.1f,%.1f)\n Center of Mass @ (x,y) = (%.1f,%.1f)\nBlue crosshairs @ centroid.\n Green crosshairs @ center of mass.\nVolume=%.2fmm^3', meanGL,sdGL,numberofPixels1,numberofPixels2c,perimeter,centroid(1),centroid(2),centerofMass(1),centerofMass(2),volume);
    msgbox(message);

elseif value>3.5 & value<4.5;
    d_imadjust = imadjust(d);
    d_histeq = histeq(d);
    d_adapthisteq = adapthisteq(d);
    imshow(d_imadjust);

    axis on;
    title('Image 4');
    set(gcf,'Position',get(0,'Screensize')); %Maximizing the figure
    message =sprintf('(Image 4)Left click the mouse to begin drawing.\n Stop holding the mouse button to finish');
    uiwait(msgbox(message));
    ROI_d = imfreehand; %Line of code from Andrew that does the drawing
    binaryImage = ROI_d.createMask; %This is creating a mask from the ROI
    xy=ROI_d.getPosition;
    subplot(2,3,1);
    imshow(d_imadjust,[]);
    axis on;
    drawnow;
    title('Image 4');
    subplot(2,3,2);
    imshow(binaryImage);
    axis on;
    title('Binary mask of the Heart image');
    labeledImage=bwlabel(binaryImage);
    measurements = regionprops(binaryImage,d_imadjust, 'area','Centroid','WeightedCentroid','Perimeter');
    area = measurements.Area;
    centroid = measurements.Centroid;
    centerofMass = measurements.WeightedCentroid;
    perimeter = measurements.Perimeter;
    numberofPixels1=sum(binaryImage(:))
    numberofPixels2d=bwarea(binaryImage)
    structBoundaries = bwboundaries(binaryImage);
    xy= structBoundaries{1}; % This gives us a n by 2 arry of x and y coordinates
    x = xy(:,2); % Column of interest
    y = xy(:,1); % row of interest
    subplot(2,3,1); % plotting over the original image
    hold on;
    plot(x,y,'LineWidth',2);
    drawnow; % forcing matlab to draw it asap
    burnedImage = d_imadjust;
    burnedImage(binaryImage) = 255;
    subplot(2,3,3);
    imshow(burnedImage);
    axis on;
    caption = sprintf('New image with mask \n burned into image');
    title(caption);
    blackMaskedImage=d_imadjust;
    blackMaskedImage(~binaryImage) = 0;
    subplot(2,3,4);
    imshow(blackMaskedImage);
    axis on;
    title('Masked Outside Region');
    meanGL=mean(blackMaskedImage(binaryImage));
    sdGL = std(double(blackMaskedImage(binaryImage)));
    hold on;
    plot(centroid(1),centroid(2),'b+','MarkerSize',20,'LineWidth',2);
    plot(centerofMass(1),centerofMass(2),'g+','MarkerSize',10,'LineWidth',2);
    insideMasked = d_imadjust;
    insideMasked(binaryImage) = 0;
    subplot(2,3,5);
    imshow(insideMasked);
    axis on;
    title('Masked inside region');
    leftColumn = min(x);
    rightColumn = max(x);
    topLine = min(y);
    bottomLine = max(y);
    width = rightColumn - leftColumn + 1;
    height = bottomLine - topLine + 1;
    croppedImage = imcrop(blackMaskedImage,[leftColumn,topLine,width,height]);
    subplot(2,3,6);
    imshow(croppedImage);
    axis 'on'
    title('Cropped Image');
    hold on;
    plot(centroid(1)-leftColumn,centroid(2)-topLine,'b+','MarkerSize',20,'LineWidth',2);
    plot(centerofMass(1)-leftColumn,centerofMass(2)-topLine,'g+','MarkerSize',20,'LineWidth',2);
    volume=(.625^2)*4*numberofPixels1;
    message=sprintf('Mean value of ROI = %.3f\n SD of ROI = %.3f\nNumber of pixels =%d\nArea in pixels=%.2f\nPerimeter = %.2f\nCentroid @ (x,y) = (%.1f,%.1f)\n Center of Mass @ (x,y) = (%.1f,%.1f)\nBlue crosshairs @ centroid.\n Green crosshairs @ center of mass.\nVolume=%.2fmm^3', meanGL,sdGL,numberofPixels1,numberofPixels2d,perimeter,centroid(1),centroid(2),centerofMass(1),centerofMass(2),volume);
    msgbox(message);
    

end

