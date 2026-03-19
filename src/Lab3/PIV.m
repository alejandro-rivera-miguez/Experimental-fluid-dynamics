clear all; clc; close all;

% Taking in the images and playing with them (Putting the origin on the bottom left corner)

%Cam 1 
%A picture

imageA_full = f_readB16(append('Cam1_0001A','.b16'));
imageA = imageA_full(540:end, 1280:end);

%B picture
imageB_full=f_readB16(append('Cam1_0001B','.b16'));
imageB = imageB_full(540:end, 1280:end);

[xmax, ymax] = size(imageA_full);

% Window sizes
wsize = [32, 32];
% wsize = [64, 64];
% wsize = [128, 128];
w_width = wsize(1);
w_height = wsize(2);

% Center points grid
xmin = w_width / 2;
ymin = w_height / 2;
xgrid = 33:w_width/2:1589; %1621-32
ygrid = 33:w_height/2:1249; %1281-32
% xgrid = 65:w_width/2:1557; %1621-64
% ygrid = 65:w_height/2:1217; %1281-64
% xgrid = 129:w_width/2:1493; %1621-128
% ygrid = 129:w_height/2:1153; %1281-128

% Number of windows in total
w_xcount = length(xgrid);
w_ycount = length(ygrid);

% These correspond to the ranges for "search" windows in image B
x_disp_max = w_width / 2;
y_disp_max = w_height / 2;


% For every window, first we have to "create" the test matrix in image A.
% Then in image B, we have to correlate this test window around its
% original position in image A, the range is pre-determined.
% The point of max. correlation corresponds to the final avg. displacement of that window.

test_ima(w_width, w_height) = 0;
test_imb(w_width + 2 * x_disp_max, w_height + 2 * y_disp_max) = 0;

dpx(w_xcount, w_ycount) = 0;
dpy(w_xcount, w_ycount) = 0;

xpeak1 = 0;
ypeak1 = 0;

% i, j are for the windows
% test_i and test_j are for the test window to be extracted from image A
for i = 1:(w_xcount)
    for j = 1:(w_ycount)
        max_correlation = 0;
        test_xmin = xgrid(i) - w_width / 2;
        test_xmax = xgrid(i) + w_width / 2;
        test_ymin = ygrid(j) - w_height / 2;
        test_ymax = ygrid(j) + w_height / 2;

        x_disp = 0;
        y_disp = 0;

        test_ima = imageA(test_xmin:test_xmax, test_ymin:test_ymax);
        test_imb = imageB((test_xmin - x_disp_max):(test_xmax + x_disp_max), ...
                          (test_ymin - y_disp_max):(test_ymax + y_disp_max));

        correlation = normxcorr2(test_ima, test_imb);
        [xpeak, ypeak] = find(correlation == max(correlation(:)));

        % Re-scaling
        xpeak1 = test_xmin + xpeak - wsize(1) / 2 - x_disp_max;
        ypeak1 = test_ymin + ypeak - wsize(2) / 2 - y_disp_max;

        dpx(i, j) = xpeak1 - xgrid(i);
        dpy(i, j) = ypeak1 - ygrid(j);
    end
end

% Vector Display
quiver(dpy, -dpx);
