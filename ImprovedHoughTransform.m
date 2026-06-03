%% Improved Hough transform
% Implementing the gradient-weighted Hough transform (GWHT), with potential 
% improvements and refinements.
% 
% By default, this program uses _testImage.png_. It can use an example ISAR 
% image instead, if _useISAR_ is set to true.

useISAR = false;
%% Load image

I = imread('testImage.png');
I = im2gray(I);
I = double(rescale(I));
figure()
imshow(I); title('Original image')
%% Image preparation
% Add some blur and noise

rng(123579); % Set fixed seed
Igauss = imgaussfilt(I,sqrt(2), 'FilterSize',5);
% Add default noise
Inoise = imnoise(Igauss);
Ispeckle = imnoise(Igauss,'speckle');
figure()
imshow(Inoise); title('Image w/ blur and noise')
% Calculating gradients

% Magnitude and direction

% If using ISAR image, skip here
if useISAR
    load("ISARexample.mat")
    [Gxa, Gya] = gradientByRatio_v2(image, 20);
else
    [Gxa, Gya] = gradientByRatio_v2(Inoise, 20);
end
Imagnitude = hypot(Gxa,Gya);
Iorientation = rad2deg(atan2(Gya,Gxa));

if useISAR
    Imagnitude(Imagnitude<0.01)=NaN; % thresholding to remove background
end

figure(); tiledlayout();
nexttile
imagesc(Imagnitude); colormap turbo; axis image; axis off; title('Magnitude')
nexttile
imagesc(Iorientation); colormap turbo; axis image; axis off; title('Orientation')
%% Thresholding

% Thresholding
imScaledGradient = rescale(Imagnitude); % Note that rescale doesn't work if the input array contains NaNs
% otsuLevel = graythresh(imScaledGradient);
otsuLevels = multithresh(imScaledGradient,2);
imOtsu = imbinarize(imScaledGradient,otsuLevels(1));
%% 
% Compare otsu levels

figure();imshowpair(imbinarize(imScaledGradient,otsuLevels(1)), ...
    imbinarize(imScaledGradient,otsuLevels(2))); title('Comparison of thresholds')

%%% Would be good to reimplement my Canny hysteresis processing here
%% 
% IMPLEMENTING CANNY HYSTERESIS HERE
% 
% BLAHBLAHBLAH


% Region reduction
CC = bwconncomp(imOtsu);
regionSizes = cellfun(@length,CC.PixelIdxList);
largeEnoughIdx = find(regionSizes>100);
p = regionprops(CC,"Area");
[~,maxIdx] = maxk([p.Area],2);
BW2 = cc2bw(CC,"ObjectsToKeep",maxIdx);

% Masks
mask = imOtsu;
maskedOrientation   = immultiply(Iorientation,mask); % Take only the data from the thresholded image
maskedMagnitude     = immultiply(Imagnitude,mask);
maskedOrientation(maskedOrientation==0) =NaN;
maskedMagnitude(maskedMagnitude==0)     =NaN;

figure(); tiledlayout(1,2); title('Masked')
nexttile
h1 = imagesc(maskedMagnitude); colormap turbo; axis image; axis off; title('Magnitude')
set(h1, 'AlphaData', ~isnan(maskedMagnitude))
nexttile
h2 = imagesc(maskedOrientation); colormap turbo; axis image; axis off; title('Orientation')
set(h2, 'AlphaData', ~isnan(maskedOrientation))
%% Hough transform

sHough = GWHT(maskedOrientation,maskedMagnitude);
%% 
% 

figure();
imagesc(sHough.hough, 'XData',sHough.theta, 'YData',sHough.rho);
colormap hsv; colorbar
xlabel('\theta'); ylabel('\rho'); 
title('GWHT')
%% Extend Hough space
% Hough space is now inherently extended to [-210,210] degrees within <./GWHT.m 
% GWHT.m>.
% 
% Peak relocating and clustering will happen after peaks are found.
%% Finding peaks

integerHough = round(sHough.hough.*100); %Multiplication to improve accuracy of rounding
if useISAR
    peakThresh = 0.05*max(integerHough(:));
else
    peakThresh = 0.1*max(integerHough(:));
end

rhoSpace = sHough.rho;
numPeaks = 100;
peaks = houghpeaks(integerHough, ...
    numPeaks, ...
    "Threshold", ceil(peakThresh));%,...
    % "NHoodSize", [51,51]);

%% 
% So a suppression neighbourhood isn't necessarily the best thing, since it 
% means that "peaks" on the edges of the cut out zone will be considered (when 
% they're not really peaks)
% 
% Instead, I could develop/implement some peak-climbing algorithm
%% Peak optimisation
% The goal here is to find initially many peaks, and then cut them down to the 
% significant ones.
% 
% Find the magnitude of each peak

linearIdx = sub2ind(size(sHough.hough), peaks(:,1), peaks(:,2));
values = sHough.hough(linearIdx)
peaks(:,3) = values;
% Sort peaks in ascending theta order
peaks = sortrows(peaks,2)
%% 
% _*peaks_ is [rho, theta, magnitude]*
%% Plotting peaks

figure();
imagesc(sHough.hough, 'XData',sHough.theta, 'YData',sHough.rho);
colormap hsv; colorbar
xlabel('\theta'); ylabel('\rho'); 
title('GWHT (w/ peaks)')
x = sHough.theta(peaks(:,2)); y = sHough.rho(peaks(:,1));
hold on
plot(x,y,'s','color','white');
hold off
%% 
% _WORK IN PROGRESS_
%% Peak clustering
% Using DBSCAN to associate detected peaks which belong to the same true peak 
% together.
% 
% Due to the extended theta domain, some peaks may be doubled. These should 
% be clustered separately, and then have clusters matched together.
% 
% Consider a weighted averaging to find the "true" peak.
% 
% Set minpts = 1, because a single detection is also valid.
% 
% How do we determine epsilon? Tricky. Consider using different distance measurements.

if useISAR
    eps = 100; % Initial test
else
    eps = 90;
end
minPts = 1;
clusterLabel = dbscan(peaks, eps, minPts);
clustColours = hsv(length(unique(clusterLabel)));
colourlist = zeros([height(peaks) 3]);
for i=1:length(colourlist)
    if clusterLabel(i) ~= -1
        colourlist(i,:) = clustColours(clusterLabel(i),:);
    else % Specific colour for ungrouped bozos
        colourlist(i,:) = [0,0,0];
    end
end
length(unique(clusterLabel))
figure();
gscatter(peaks(:,2), peaks(:,1), clusterLabel)
%% 
% Eps should be refined for each purpose - ideally automated!
%% Associate clusters
% Thanks to the expanded Hough domain, paired clusters will be 360° apart

% Find theta resolution, and rho resolution (for later)
thetaRes = sHough.theta(2)-sHough.theta(1);
rhoRes = sHough.rho(2)-sHough.rho(1);
% Find how many theta vals later that corresponds to
indexDiff = 360/thetaRes;
%% 
% We only need to look at the first 60 degrees

% Similarly calculated
stopIdx = (60/thetaRes)+1; % The index at which to stop checking
%% 
% Iterate through all, until 60° has passed.

% Initialise merging pairs list
mergingPairs = [];
for i=2%1:height(sortedpeaks)
    if peaks(i,2)<stopIdx
        corrTheta = round(peaks(i,2) + indexDiff); % This has to be an integer
        corrPoint = [peaks(i,1), corrTheta];

        % Check if this exists
        idx = find(ismember(peaks(:,1:2), corrPoint, "rows"));
        if ~isempty(idx)
            % Correlate these points
            mergingPairs = [mergingPairs; clusterLabel(i), clusterLabel(idx)]
        end

    else
        break
    end
end
%% 
% Correlate the clusters

% Make a graph
if ~isempty(mergingPairs)
    G = graph(mergingPairs(:,1), mergingPairs(:,2))
    % Find connected components
    componentID = conncomp(G)
    % Initialise new label list
    connClustLabel = zeros(size(clusterLabel));
    % Now relabel points
    for i=1:length(clusterLabel)
        oldLabel = clusterLabel(i);
        newLabel = componentID(oldLabel);
        connClustLabel(i)=newLabel;
    end
else
    disp('No pairs to merge')
    connClustLabel = clusterLabel;
end
%% Finding average peak location
% For each cluster, there must be a "true" peak location. We can take a weighted 
% average.
% 
% _There will be an edge case for clusters that wrap around the boundaries. 
% Computing circular averages must be considered - there's a circular statistics 
% toolbox for MATLAB that could be good to look into._
% 
% In the general case

% For each unique cluster
meanPeaks = [];
for i=1:length(unique(connClustLabel))
    % Slice all the peaks belonging to this cluster
    clusterPeaks = peaks(connClustLabel==i,:);
    
    % Some check here to see if the cluster is disparate
    %%% Since clusterPeaks is already sorted according to theta, check the
    %%% distance between first theta and last theta
    if clusterPeaks(end,2)-clusterPeaks(1,2) > 0.5*indexDiff
        % Threshold chosen to be larger than any legitimate nondisparate cluster, but
        % smaller than any legitimate disparate cluster
        disp('Disparity detected')
        % If disparate, add to the smaller values - any theta values less
        % than 0° get 360° added to them
        for j=1:height(clusterPeaks)
            if clusterPeaks(j,2) < 210/thetaRes
                clusterPeaks(j,2) = clusterPeaks(j,2) + 360/thetaRes;
            end
        end
    end

    % Calculate a weighted mean - but weight STRONGLY towards better
    % magnitude, by rescaling and squaring.
    if height(clusterPeaks)~=1 % For clusters with only one peak, the rescaling was setting the value to 0
        weights = (rescale(clusterPeaks(:,3))).^2;
    else
        weights = 1;
    end
    meanRho     = mean(clusterPeaks(:,1),'Weights',weights);
    meanTheta   = mean(clusterPeaks(:,2),'Weights',weights);

    meanPeaks = [meanPeaks; meanRho, meanTheta];
end
%% 
% Let's plot these peaks and see how well they work on this data.

figure();
imagesc(sHough.hough, 'XData',sHough.theta, 'YData',sHough.rho);
colormap hsv; colorbar
xlabel('\theta'); ylabel('\rho'); 
title('GWHT (w/ peaks)')

x = min(sHough.theta) + meanPeaks(:,2)*thetaRes;
y = min(sHough.rho) + meanPeaks(:,1)*rhoRes;
% x = sHough.theta(peaks(:,2)); y = sHough.rho(peaks(:,1));
hold on
plot(x,y,'s','color','white');
hold off
xline(-180)
%% Reconstruct lines

for j=1:length(meanPeaks)
    pk = round(meanPeaks(j,:));
    % Find a new mask for that peak. Within a range of the peak angle
    pkAngle = sHough.theta(pk(2));
    pkMask = zeros(size(mask));
    pkMask(maskedOrientation>(pkAngle-10) & maskedOrientation<(pkAngle+10)) = 1;
    pkLine = houghlines(pkMask, sHough.theta,rhoSpace,pk,...
        "FillGap",50, "MinLength", 100);

    if j==1 % For the first one, initialise the new big struct
        sLines = pkLine;
    else
        sLines = [sLines, pkLine];
    end
end
%%
tLines = struct2table(sLines);
thetaRhoVals = [sLines.theta;sLines.rho]';
uniqueVals = unique(thetaRhoVals,"rows");
colours = prism(length(uniqueVals));

for i=1:length(uniqueVals) %for each unique value
    idx = and(any(thetaRhoVals == uniqueVals(i,1),2), ...
        any(thetaRhoVals == uniqueVals(i,2),2)); % Matching rows
    tLines.label(idx) = i;
end

figure();
if useISAR
    imshow(image,[])
else
    imshow(Inoise)
end
hold on
for k = 1:height(tLines)
   xy = [tLines.point1(k,:); tLines.point2(k,:)];
   colour = colours(tLines.label(k),:);
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color',colour);
end
title('Line detection overlay');

hold off