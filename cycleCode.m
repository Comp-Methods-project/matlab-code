clear
clc

[files,folder] = uigetfile('*.dcm','MultiSelect','on');
% Allows user to select one, or multiple dcm files
[k long]= size(files);
% Assigning value of number of files selected to variable "long"




app.Slider_2.Limits = [0 long];
app.Slider_2.MajorTicks=[0:1:long];
app.Slider_2.MinorTicks=[0 long];
if ~iscell(files)
    files={files};
end

for ii=1:length(files)
    infos(ii)=dicominfo([folder,filesep,files{ii}]);
    % infos has all the filenames, one on each row
    data(:,:,ii)=dicomread([folder,filesep,files{ii}]);
    % data is image width x height x number of files
    slider(:,:,ii)=[0 ii];
end

app.allinfos=infos;
app.alldata=data;
app.slider=slider;
            
            

count = 1;
% initializing count for the while loop below          

      
      
while count <=long   
    
    I=dicomread(fullfile(files{count}));
    % fullfile allows us to access actual filename. This is necessary for
    % dicomread command
    a_imadjust=imadjust(I);
    % adjusting image intensity of dcm files
    imshow(a_imadjust);
    % showing dcm images on the screen
    
    axis on;
    title('Original Image'); % Image title
    set(gcf,'Position',get(0,'Screensize')); % Maximizing the figure
   
    %Now asking user to draw on the image, after which we will apply a mask
    message =sprintf('Left click the mouse to begin drawing.\n Stop holding the mouse button to finish');
    uiwait(msgbox(message));
    ROI_a = imfreehand; %Line of code that allows user to select ROI
    binaryImage = ROI_a.createMask;% This is creating a mask from the ROI
    % Mask is basically a binary image of the ROI
    xy=ROI_a.getPosition;
    % After the user drew ROI, the masking line of code gave a pixel value of 1
    % to the image pixels belonging to ROI.  The rest of the pixels recieved a
    % value of 0, setting them up as part of the background
    
    subplot(2,3,1);
    imshow(a_imadjust,[]);
    axis on;
    drawnow;
    title('Original Image');
    
    % Now we display the mask that was drawn by user
    subplot(2,3,2);
    imshow(binaryImage);
    axis on;
    title('Binary mask of the image');
    
    
    nextLoop = 1;
    structBoundaries = bwboundaries(binaryImage);
    if isempty(structBoundaries);
        nextLoop = 0;
        message2=sprintf('You did not select sufficiently large ROI\nThis program will now terminate\nNo total volume calculation will be provided');
        uiwait(msgbox(message2));
        close all
        break
        
    end
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
    
    numberofPixels1=sum(binaryImage(:));
    numberofPixels2(count)=sum(binaryImage(:));
    Volume = (.625^2)*4*numberofPixels1;
    TotalVolume = (.625^2)*4*numberofPixels2;
    message=sprintf('Number of pixels =%d\nVolume of ROI=%.2f',numberofPixels1,Volume);
    uiwait(msgbox(message));
    close all
    count = count+1;
    
end


nextCount = 2;


while nextCount > nextLoop;
    if nextLoop == 0;
       nextCount = 0;
    break
else
    SumVolume = sum(TotalVolume);
    message1=sprintf('Total volume from all selected slices = %.2fmm^3',SumVolume);
    uiwait(msgbox(message1));
    end
    nextCount = nextCount - 1;
end

