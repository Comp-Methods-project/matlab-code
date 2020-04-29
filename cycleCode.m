 clear
 clc

[files,folder] = uigetfile('*.dcm','MultiSelect','on');
            long = length(files);
            if long == 16;
                long = 1;
            else
            long == long;
             end
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
                % slider is generating values that I'm not sure how to use
            end
            app.allinfos=infos;
            app.alldata=data;
            app.slider=slider;
            
            
            
         
            
  count = 1;           
            
while count <=long;
    for i = 1<=long;
    I=dicomread(fullfile(files{i}));
    a_imadjust=imadjust(I);
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
    

    numberofPixels1=sum(binaryImage(:))
    
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
Volume = (.625^2)*4*numberofPixels1;
    
    message=sprintf('Number of pixels =%d\nVolume of ROI=%.2f',numberofPixels1,Volume);
    uiwait(msgbox(message));
    close all
    end
    count = count +1
end