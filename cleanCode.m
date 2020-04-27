clear
clc

a = dicomread('IM-0001-0001.dcm'); % read DICOM file
a_imadjust = imadjust(a); % adjust the image brightness
figure(5)
imshow(a_imadjust) % show the image
axis on;
title('Original MRI image'); % Image title
set(gcf,'Position',get(0,'Screensize')); %Maximizing the figure


%Now asking user to draw on the image, after which we will apply a mask
message =sprintf('Left click the mouse to begin drawing.\n Stop holding the mouse button to finish');
uiwait(msgbox(message)); 
ROI_a = imfreehand; %Line of code that allows user to select ROI
binaryImage = ROI_a.createMask; %This is creating a mask from the ROI
%Mask is basically a binary image of the ROI
xy=ROI_a.getPosition;
% After the user drew ROI, the masking line of code gave a pixel value of 1
% to the image pixels belonging to ROI.  The rest of the pixels recieved a
% value of 0, setting them up as part of the background


%using subplots to show more images
subplot(2,3,1);
imshow(a_imadjust,[]);
axis on;
drawnow;
title('Original MRI image');


% Now we display the mask that was drawn by user
subplot(2,3,2);
imshow(binaryImage);
axis on;
title('Binary mask of the MRI image');

info1 = dicominfo('IM-0001-0001.dcm');
% This is to look at metadata of DICOM file

info1.PixelSpacing; % output in mm
% ans =
% 
%     0.6250
%     0.6250

% The first value of pixel spacing is telling us the vertical spacing
% between two adjacent pixels. 
% The second value of pixel spacing is telling us horizontal spacing
% between two adjacent pixels.

info1.SliceThickness; % output in mm
% ans =
% 
%      4

% As the name suggests, this is giving us the z value of the slice

info1.SpacingBetweenSlices; %output in mm
% ans =
% 
%      4

% Because this value is same as slice thickness, it means that there is no
% gap between any of the four slices, assuming that they are concurrent.


% Calculating the number of pixels and volume from pixels selected by user
numberofPixels1=sum(binaryImage(:));
Volume = (.625^2)*4*numberofPixels1;
% In calculating the ROI volume, we have to first calculate volume of
% each pixel. Our length = .625 mm, width = .625mm, and height = 4mm. The
% value of this calculation can be multiplied with total number of pixels
% from ROI, which gives us our volume. 


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


% Code below is  setting the pixel values of ROI to 0 in the
% original image. 
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


% report the results of the calculation
message=sprintf('Number of pixels =%d\nVolume of ROI=%.2f',numberofPixels1,Volume);
msgbox(message);