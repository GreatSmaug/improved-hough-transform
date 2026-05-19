%% Standard Hough transform implementation, to contrast with ImprovedHoughTransform

%% Load data
% Toggle useISAR to use ISAR image instead
useISAR = false;
% Load image
I = imread('testImage.png');
I = im2gray(I);
I = double(rescale(I));
figure()
imshow(I); title('Original image')

% Image preparation
% Add some blur and noise
rng(123579); % Set fixed seed
Igauss = imgaussfilt(I,sqrt(2), 'FilterSize',5);
% Add default noise
Inoise = imnoise(Igauss);
% Ispeckle = imnoise(Igauss,'speckle');
imshow(Inoise); title('Image w/ blur and noise')

%% Hough processing
% Note that the standard Hough implementation doesn't use gradient-by-ratio
% calculation. Instead we use the inbuilt edge function
BW = edge(Inoise,"canny");
imshow(BW); title('Edge image')
% Match Hough parameters as closely as possible
[H,T,R] = hough(BW,Theta=-90:0.2:89.8);

%% Display Hough space
figure();
imagesc(H, 'XData', T, 'YData', R);
colormap hsv; colorbar
xlabel('\theta'); ylabel('\rho'); 
title('SHT')

%% Find peaks
peakThresh = 0.1*max(H(:));
numPeaks = 100;
peaks = houghpeaks(H, ...
    numPeaks, ...
    "Threshold", ceil(peakThresh));
%% Plot peaks
x = T(peaks(:,2)); y = R(peaks(:,1));
hold on
plot(x,y,'s','color','white');
hold off

%% Find lines
% Again, no refined implementation
lines = houghlines(BW,T,R,peaks,"FillGap",5,"MinLength",7);

%% Plot lines
imshow(Inoise)
hold on
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),LineWidth=2,Color="green");
   plot(xy(1,1),xy(1,2),"x",LineWidth=2,Color="yellow");
   plot(xy(2,1),xy(2,2),"x",LineWidth=2,Color="red");
end
title('Line detection overlay');