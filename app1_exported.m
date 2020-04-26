classdef app1_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                  matlab.ui.Figure
        ViewSeparateImagesButton  matlab.ui.control.Button
        ViewStackedImagesButton   matlab.ui.control.Button
        ImageSelectorSliderLabel  matlab.ui.control.Label
        ImageSelectorSlider       matlab.ui.control.Slider
        Image                     matlab.ui.control.Image
        Image1Label               matlab.ui.control.Label
        Image2                    matlab.ui.control.Image
        Image2Label               matlab.ui.control.Label
        Image3                    matlab.ui.control.Image
        Image3Label               matlab.ui.control.Label
        Image4                    matlab.ui.control.Image
        Image4Label               matlab.ui.control.Label
        ToselectacertainImageusesliderbelowLabel  matlab.ui.control.Label
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: ViewSeparateImagesButton
        function ViewSeparateImagesButtonPushed(app, event)
a = dicomread('IM-0001-0001.dcm');
b = dicomread('IM-0001-0002.dcm');
c = dicomread('IM-0001-0003.dcm');
d = dicomread('IM-0001-0004.dcm');

figure(1)
a_imadjust = imadjust(a);
a_histeq = histeq(a);
a_adapthisteq = adapthisteq(a);
imshow(a_imadjust)

figure(2)
b_imadjust = imadjust(b);
b_histeq = histeq(b);
b_adapthisteq = adapthisteq(b);
imshow(b_imadjust)

figure(3)
c_imadjust = imadjust(c);
c_histeq = histeq(c);
c_adapthisteq = adapthisteq(c);
imshow(c_imadjust)

figure(4)
d_imadjust = imadjust(d);
d_histeq = histeq(d);
d_adapthisteq = adapthisteq(d);
imshow(d_imadjust)
        end

        % Button pushed function: ViewStackedImagesButton
        function ViewStackedImagesButtonPushed(app, event)
a = dicomread('IM-0001-0001.dcm');
b = dicomread('IM-0001-0002.dcm');
c = dicomread('IM-0001-0003.dcm');
d = dicomread('IM-0001-0004.dcm');

D = cat(3,a,b,c,d);
figure
colormap gray
contourslice(D,[],[],[1,2,3,4,27],15);
view(3)
axis tight

        end

        % Value changed function: ImageSelectorSlider
        function ImageSelectorSliderValueChanged(app, event)
            value = app.ImageSelectorSlider.Value;
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
    message=sprintf('Mean value of ROI = %.3f\n SD of ROI = %.3f\nNumber of pixels =%d\nArea in pixels=%.2f\nPerimeter = %.2f\nCentroid @ (x,y) = (%.1f,%.1f)\n Center of Mass @ (x,y) = (%.1f,%.1f)\nBlue crosshairs @ centroid.\n Green crosshairs @ center of mass.\nVolume=%.2fmm^3', meanGL,sdGL,numberofPixels1,numberofPixels2d,perimeter,centroid(1),centroid(2),centerofMass(1),centerofMass(2),volume);
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
    volume=(.625^2)*4*numberofPixels1;
    message=sprintf('Mean value of ROI = %.3f\n SD of ROI = %.3f\nNumber of pixels =%d\nArea in pixels=%.2f\nPerimeter = %.2f\nCentroid @ (x,y) = (%.1f,%.1f)\n Center of Mass @ (x,y) = (%.1f,%.1f)\nBlue crosshairs @ centroid.\n Green crosshairs @ center of mass.\nVolume=%.2fmm^3', meanGL,sdGL,numberofPixels1,numberofPixels2d,perimeter,centroid(1),centroid(2),centerofMass(1),centerofMass(2),volume);
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
    message=sprintf('Mean value of ROI = %.3f\n SD of ROI = %.3f\nNumber of pixels =%d\nArea in pixels=%.2f\nPerimeter = %.2f\nCentroid @ (x,y) = (%.1f,%.1f)\n Center of Mass @ (x,y) = (%.1f,%.1f)\nBlue crosshairs @ centroid.\n Green crosshairs @ center of mass.\nVolume=%.2fmm^3', meanGL,sdGL,numberofPixels1,numberofPixels2d,perimeter,centroid(1),centroid(2),centerofMass(1),centerofMass(2),volume);
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
    volume=(.625^2)*4*numberofPixels1;
    message=sprintf('Mean value of ROI = %.3f\n SD of ROI = %.3f\nNumber of pixels =%d\nArea in pixels=%.2f\nPerimeter = %.2f\nCentroid @ (x,y) = (%.1f,%.1f)\n Center of Mass @ (x,y) = (%.1f,%.1f)\nBlue crosshairs @ centroid.\n Green crosshairs @ center of mass.\nVolume=%.2fmm^3', meanGL,sdGL,numberofPixels1,numberofPixels2d,perimeter,centroid(1),centroid(2),centerofMass(1),centerofMass(2),volume);
    msgbox(message);
    

end



        end

        % Image clicked function: Image
        function ImageClicked(app, event)
a = dicomread('IM-0001-0001.dcm');
a_imadjust = imadjust(a);
a_histeq = histeq(a);
a_adapthisteq = adapthisteq(a);
imshow(a_imadjust)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'UI Figure';

            % Create ViewSeparateImagesButton
            app.ViewSeparateImagesButton = uibutton(app.UIFigure, 'push');
            app.ViewSeparateImagesButton.ButtonPushedFcn = createCallbackFcn(app, @ViewSeparateImagesButtonPushed, true);
            app.ViewSeparateImagesButton.Position = [61 380 135 22];
            app.ViewSeparateImagesButton.Text = 'View Separate Images';

            % Create ViewStackedImagesButton
            app.ViewStackedImagesButton = uibutton(app.UIFigure, 'push');
            app.ViewStackedImagesButton.ButtonPushedFcn = createCallbackFcn(app, @ViewStackedImagesButtonPushed, true);
            app.ViewStackedImagesButton.Position = [277 380 131 22];
            app.ViewStackedImagesButton.Text = 'View Stacked Images';

            % Create ImageSelectorSliderLabel
            app.ImageSelectorSliderLabel = uilabel(app.UIFigure);
            app.ImageSelectorSliderLabel.HorizontalAlignment = 'right';
            app.ImageSelectorSliderLabel.FontWeight = 'bold';
            app.ImageSelectorSliderLabel.Position = [86 87 92 22];
            app.ImageSelectorSliderLabel.Text = 'Image Selector';

            % Create ImageSelectorSlider
            app.ImageSelectorSlider = uislider(app.UIFigure);
            app.ImageSelectorSlider.Limits = [0 4];
            app.ImageSelectorSlider.MajorTicks = [0 1 2 3 4];
            app.ImageSelectorSlider.ValueChangedFcn = createCallbackFcn(app, @ImageSelectorSliderValueChanged, true);
            app.ImageSelectorSlider.MinorTicks = [1 2 3 4];
            app.ImageSelectorSlider.FontWeight = 'bold';
            app.ImageSelectorSlider.Position = [199 96 150 3];

            % Create Image
            app.Image = uiimage(app.UIFigure);
            app.Image.ImageClickedFcn = createCallbackFcn(app, @ImageClicked, true);
            app.Image.BackgroundColor = [0.651 0.651 0.651];
            app.Image.Position = [32 214 100 100];
            app.Image.ImageSource = '1.png';

            % Create Image1Label
            app.Image1Label = uilabel(app.UIFigure);
            app.Image1Label.BackgroundColor = [0 1 1];
            app.Image1Label.FontWeight = 'bold';
            app.Image1Label.Position = [57 333 54 22];
            app.Image1Label.Text = ' Image 1';

            % Create Image2
            app.Image2 = uiimage(app.UIFigure);
            app.Image2.BackgroundColor = [0.651 0.651 0.651];
            app.Image2.Position = [131 214 100 100];
            app.Image2.ImageSource = '2.png';

            % Create Image2Label
            app.Image2Label = uilabel(app.UIFigure);
            app.Image2Label.BackgroundColor = [0 1 1];
            app.Image2Label.FontWeight = 'bold';
            app.Image2Label.Position = [156 333 54 22];
            app.Image2Label.Text = ' Image 2';

            % Create Image3
            app.Image3 = uiimage(app.UIFigure);
            app.Image3.BackgroundColor = [0.651 0.651 0.651];
            app.Image3.Position = [230 214 100 100];
            app.Image3.ImageSource = '3.png';

            % Create Image3Label
            app.Image3Label = uilabel(app.UIFigure);
            app.Image3Label.BackgroundColor = [0 1 1];
            app.Image3Label.FontWeight = 'bold';
            app.Image3Label.Position = [255 333 54 22];
            app.Image3Label.Text = ' Image 3';

            % Create Image4
            app.Image4 = uiimage(app.UIFigure);
            app.Image4.BackgroundColor = [0.651 0.651 0.651];
            app.Image4.Position = [329 214 100 100];
            app.Image4.ImageSource = '4.png';

            % Create Image4Label
            app.Image4Label = uilabel(app.UIFigure);
            app.Image4Label.BackgroundColor = [0 1 1];
            app.Image4Label.FontWeight = 'bold';
            app.Image4Label.Position = [354 333 54 22];
            app.Image4Label.Text = ' Image 4';

            % Create ToselectacertainImageusesliderbelowLabel
            app.ToselectacertainImageusesliderbelowLabel = uilabel(app.UIFigure);
            app.ToselectacertainImageusesliderbelowLabel.BackgroundColor = [0 1 1];
            app.ToselectacertainImageusesliderbelowLabel.FontWeight = 'bold';
            app.ToselectacertainImageusesliderbelowLabel.Position = [131 167 262 22];
            app.ToselectacertainImageusesliderbelowLabel.Text = ' To select a certain Image, use slider below: ';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = app1_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

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