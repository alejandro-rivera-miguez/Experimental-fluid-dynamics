%%
clear all, close all, clc
%% Image format convert. Converts all of the .b16 images given in Lab3 into
% .png format for treatment with PIVLab. Makes use of the ready made
% f_readB16.m function and the image processing toolbox. Basic workflow
% tested in Format_Conversion_Trial.m

%Essentially a for loop acting on a series of files. Note there are files
%for cameras 1 and 2 and for pictures A and B in each of the 100 PIV
%frames. The for loop itself is divided to avoid multiple if statements to
%manage the format selection from the given files. 

for i=1:9
    %Cam 1 
    %A picture
    str=sprintf('Gruppo5/10ms/Cam1/Cam1_000%dA',i);
    image=f_readB16(append(str,'.b16'));
    image=uint16(image);
    imwrite(image,append(str,'.png'));
    %B picture
    str=sprintf('Gruppo5/10ms/Cam1/Cam1_000%dB',i);
    image=f_readB16(append(str,'.b16'));
    image=uint16(image);
    imwrite(image,append(str,'.png'));
    %Cam 2 
    %A picture
    str=sprintf('Gruppo5/10ms/Cam2/Cam2_000%dA',i);
    image=f_readB16(append(str,'.b16'));
    image=uint16(image);
    imwrite(image,append(str,'.png'));
    %B picture
    str=sprintf('Gruppo5/10ms/Cam2/Cam2_000%dB',i);
    image=f_readB16(append(str,'.b16'));
    image=uint16(image);
    imwrite(image,append(str,'.png'));
end

for i=10:99
    %Cam 1 
    %A picture
    str=sprintf('Gruppo5/10ms/Cam1/Cam1_00%dA',i);
    image=f_readB16(append(str,'.b16'));
    image=uint16(image);
    imwrite(image,append(str,'.png'));
    %B picture
    str=sprintf('Gruppo5/10ms/Cam1/Cam1_00%dB',i);
    image=f_readB16(append(str,'.b16'));
    image=uint16(image);
    imwrite(image,append(str,'.png'));
    %Cam 2 
    %A picture
    str=sprintf('Gruppo5/10ms/Cam2/Cam2_00%dA',i);
    image=f_readB16(append(str,'.b16'));
    image=uint16(image);
    imwrite(image,append(str,'.png'));
    %B picture
    str=sprintf('Gruppo5/10ms/Cam2/Cam2_00%dB',i);
    image=f_readB16(append(str,'.b16'));
    image=uint16(image);
    imwrite(image,append(str,'.png'));
end

i=100;
%Cam 1 
%A picture
str=sprintf('Gruppo5/10ms/Cam1/Cam1_0%dA',i);
image=f_readB16(append(str,'.b16'));
image=uint16(image);
imwrite(image,append(str,'.png'));
%B picture
str=sprintf('Gruppo5/10ms/Cam1/Cam1_0%dB',i);
image=f_readB16(append(str,'.b16'));
image=uint16(image);
imwrite(image,append(str,'.png'));
%Cam 2 
%A picture
str=sprintf('Gruppo5/10ms/Cam2/Cam2_0%dA',i);
image=f_readB16(append(str,'.b16'));
image=uint16(image);
imwrite(image,append(str,'.png'));
%B picture
str=sprintf('Gruppo5/10ms/Cam2/Cam2_0%dB',i);
image=f_readB16(append(str,'.b16'));
image=uint16(image);
imwrite(image,append(str,'.png'));