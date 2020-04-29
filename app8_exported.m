classdef app8_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                matlab.ui.Figure
        SelectROIVOIButton      matlab.ui.control.Button
        PushforhappinessButton  matlab.ui.control.Button
        AblatedAtrialTissueVolumeCalculatorLabel  matlab.ui.control.Label
        Image                   matlab.ui.control.Image
    end

    
    properties (Access = private)
        alldata;
        allinfos;% Description
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, a)
clear all
clc
close all
        end

        % Button pushed function: SelectROIVOIButton
        function SelectROIVOIButtonPushed(app, event)
clear
clc
close all

% Allows user to select one, or multiple dcm files
files = cellstr(uigetfile('*.dcm','MultiSelect','on'));
[k long]= size(files);
% Assigning value of number of columns from files cell array to variable
% "long"

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
    
    info=dicominfo(fullfile(files{count}));
    % extracting out metadata from DCM file
    
    PixelSize=info.PixelSpacing;
    % retrieving Pixel Spacing info from image file
    
    Pixel_length = PixelSize(1);
    % Assigning pixel length info to variable "Pixel_length"; in mm.
    
    Pixel_width = PixelSize(2);
    % Assigning pixel width info to variable "Pixel_width"; in mm.
    
    PixelThickness = info.SliceThickness;
    % Retrieving info about thickness of the slice from DCM metadata. This
    % info will function as z axis value for pixel. Output is in mm.   
    
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
    % Variable initialized to keep count for "total volume calculation
    % while loop"
    
    structBoundaries = bwboundaries(binaryImage);
    % Function of "if" conditional below is twofold
    % 1) If the user fails to select a ROI by just clicking and not drawing
    % any ROI, the conditional below will initialize and halt the entire
    % program. 
    % 2) In order to make sure that the total volume while loop doesn't
    % error out, this conditional will assign a value of 0 to
    % aforementioned while loop, ensuring it does not execute in the first
    % place
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
    
    % Now keeping only the part of the image that is inside the mask
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
    
    % Counting the number of pixels. We are storing the count in two
    % separate variables. The first variable will have its counting
    % displayed for each slice output. The second variable is going to keep
    % track of number of pixels from all the slices selected in a cell
    % array
    numberofPixels1=sum(binaryImage(:));
    numberofPixels2(count)=sum(binaryImage(:));
   
    % Two separate volume calculations. First volume calculation is unique
    % to each slice selected by user. The second volume calculation is
    % keeping track of all the slice volumes generated during the loop.
    % Volume calculation for both is done by calculating pixel volume, info
    % for which was retrieved drom DCM metadata earlier, and multiplying
    % that value by number of pixels from user defined ROI. 
    Volume = Pixel_length*Pixel_width*PixelThickness*numberofPixels1;
    TotalVolume = Pixel_length*Pixel_width*PixelThickness*numberofPixels2;
    
    % Number of pixels and volume output for each selected slice
    message=sprintf('Number of pixels =%d\nVolume of ROI=%.2fmm^3',numberofPixels1,Volume);
    uiwait(msgbox(message));
    close all
    
    %Moving onto next selected slice
    count = count+1;
    
end

% Counter for while loop below
nextCount = 2;


while nextCount > nextLoop;
    % "If" conditional below keeps track whether or not the program was
    % halted in the while loop above. If the while loop was halted, then
    % this while loop will halt as well, preventing calculation of total
    % volume, and its output to the user
    if nextLoop == 0;
       nextCount = 0;
    break
else
    SumVolume = sum(TotalVolume);
    message1=sprintf('Total volume from all selected slices = %.2fmm^3',SumVolume);
    uiwait(msgbox(message1));
    end
    
    % Stopping this while loop
    nextCount = nextCount - 1;
end
        end

        % Button pushed function: PushforhappinessButton
        function PushforhappinessButtonPushed(app, event)
clc
close all
clear all
x1=[1 cos(pi/12) cos(pi/6) cos(pi/4) cos(pi/3) cos(5*pi/12) cos(pi/2) -cos(5*pi/12) -cos(pi/3) -cos(pi/4) -cos(pi/6) -cos(pi/12) -1];
l=length(x1);
for i=1:l
    y1(i)=sqrt(1-(x1(i)^2));
end
plot(x1,y1)
yr1=-y1;
hold on
plot(x1,yr1)
% xx=0:0.01:1;
% yy = spline(x,y,xx);
% plot(x,y,'o',xx,yy);
% % plot(x,y,xx)
x2=0.02*x1-0.4;
y2=0.02*y1+0.4;
yr2=-y2+0.8;
plot(x2,y2)
plot(x2,yr2)
x3=0.02*x1+0.4;
y3=0.02*y1+0.4;
yr3=-y3+0.8;
plot(x3,y3)
plot(x3,yr3)
x4=[0 0.05 0.1 0.15 0.2 0.25 0.3];
l1=length(x4);
for i=1:l1
    y4(i)=2*(x4(i)^2)-0.4;
end
plot(x4,y4)
x5=-x4;
plot(x5,y4)
x6=[0 0 0];
y6=[-0.1 0 0.1];
plot(x6,y6)
axis([-2 2 -2 2]);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 381 320];
            app.UIFigure.Name = 'UI Figure';

            % Create SelectROIVOIButton
            app.SelectROIVOIButton = uibutton(app.UIFigure, 'push');
            app.SelectROIVOIButton.ButtonPushedFcn = createCallbackFcn(app, @SelectROIVOIButtonPushed, true);
            app.SelectROIVOIButton.Position = [142 74 100 22];
            app.SelectROIVOIButton.Text = 'Select ROI/VOI';

            % Create PushforhappinessButton
            app.PushforhappinessButton = uibutton(app.UIFigure, 'push');
            app.PushforhappinessButton.ButtonPushedFcn = createCallbackFcn(app, @PushforhappinessButtonPushed, true);
            app.PushforhappinessButton.Position = [133 31 119 22];
            app.PushforhappinessButton.Text = 'Push for happiness';

            % Create AblatedAtrialTissueVolumeCalculatorLabel
            app.AblatedAtrialTissueVolumeCalculatorLabel = uilabel(app.UIFigure);
            app.AblatedAtrialTissueVolumeCalculatorLabel.FontSize = 16;
            app.AblatedAtrialTissueVolumeCalculatorLabel.FontWeight = 'bold';
            app.AblatedAtrialTissueVolumeCalculatorLabel.Position = [36 272 312 22];
            app.AblatedAtrialTissueVolumeCalculatorLabel.Text = 'Ablated Atrial Tissue Volume Calculator';

            % Create Image
            app.Image = uiimage(app.UIFigure);
            app.Image.Position = [100 114 184 134];
            app.Image.ImageSource = 'ear_ecg.jpg';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = app8_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end