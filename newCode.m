clear
clc

% dicomreadVolume is reading dicom volume from input source. 
% uigetfile is allowing us to select all the dcm files in the current
% folder. Multiselect is on. Thus, mutliple files can be selected.

% dicomreadVolume only works with 2 or more slices
[V,sp,dim] = dicomreadVolume(uigetfile('*.dcm','MultiSelect','on'));

% Example V, in which 3 slices are selected
%  V = 368x512x1x3            
%  V = [rows, columns, channels, slices]
% Channels is number of color channels per voxel. Because our images are
% grayscale, our channel value = 1.

% sp = 
% 
%   struct with fields:
% 
%        PatientPositions: [3×3 double]
%           PixelSpacings: [3×2 double]
%     PatientOrientations: [2×3×3 double]


% sp contains spatial information about the voxels

% dim = dimension. With three slices, dimension would be 3. 

V = squeeze(V);
% This is getting rid of the third dimension, also known as singleton
% dimension. Voxels do not have multiple color components. 
% Now, assuming 3 slices are selected, V = 368x512x3

count = 0;


    



%% keep count
clear
clc 

filename=uigetfile('*.dcm','MultiSelect','on');
files_length=length(filename);

if files_length == 16;
    files_length = 1;
else
    files_length == files_length;
end

count = 1;


while count <= files_length;
    for i=1:files_length;
        file=dicominfo(fullfile(filename{i}));
        prop_PixelSpacing = file.PixelSpacing;
        prop_PatientOrientation = file.ImageOrientationPatient
        prop_PatientPosition = file.ImagePositionPatient
        
    end
    count = count+1;
end