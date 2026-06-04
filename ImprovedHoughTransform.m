%% Improved Hough transform
% Implementing the gradient-weighted Hough transform (GWHT), with potential 
% improvements and refinements.
% 
% By default, this program uses _testImage.png_. It can use an example ISAR 
% image instead, if _useISAR_ is set to true.

useISAR = true;
%% Load image

I = imread('testImage.png');
I = im2gray(I);
I = double(rescale(I));
figure()
imshow(I); title('Original image')
%% 
% 
% 
% Scaling command: f=gcf; f.Position = [281.8000 252.2000 424.8000 416];
%% Image preparation
% Add some blur and noise

rng(123579); % Set fixed seed
Igauss = imgaussfilt(I,sqrt(2), 'FilterSize',5);
% Add default noise
Inoise = imnoise(Igauss);
Ispeckle = imnoise(Igauss,'speckle');
fNoise = figure('Position',[281.8000 252.2000 424.8000 416]);
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
imThresh = multithresh(imScaledGradient,2)
imOtsu = imbinarize(imScaledGradient,imThresh(1));
%% 
% Compare otsu levels

figure();imshowpair(imbinarize(imScaledGradient,imThresh(1)), ...
    imbinarize(imScaledGradient,imThresh(2))); title('Comparison of thresholds')
%% Canny hysteresis

strongEdge = imScaledGradient>imThresh(2);
weakEdge = imScaledGradient>imThresh(1);
[rstrong, cstrong] = find(strongEdge); % Row and column of strong pixels
hysteresisMask = bwselect(weakEdge, cstrong, rstrong, 8);
cannyMaskedImage = immultiply(imScaledGradient,hysteresisMask);
figure();
h1 = imagesc(cannyMaskedImage); colormap turbo; axis image; axis off; 
title('Edges after hysteresis')
set(h1, 'AlphaData', ~isnan(cannyMaskedImage))
c=colorbar;c.Label.String = 'Gradient magnitude';

%% 
% 

% % Region reduction
% CC = bwconncomp(hysteresisMask);
% regionSizes = cellfun(@length,CC.PixelIdxList);
% largeEnoughIdx = find(regionSizes>100);
% p = regionprops(CC,"Area");
% [~,maxIdx] = maxk([p.Area],2);
% BW2 = cc2bw(CC,"ObjectsToKeep",maxIdx);
%% 
% Region reduction no longer needed after hysteresis

% Masks
mask = hysteresisMask;
maskedOrientation   = immultiply(Iorientation,mask); % Take only the data from the thresholded image
maskedMagnitude     = immultiply(Imagnitude,mask);
maskedOrientation(maskedOrientation==0) =NaN;
maskedMagnitude(maskedMagnitude==0)     =NaN;

figure(); tiledlayout("TileSpacing","tight", "Padding","tight");
nexttile
h1 = imagesc(maskedMagnitude); colormap turbo; axis image; axis off; title('Magnitude')
set(h1, 'AlphaData', ~isnan(maskedMagnitude))
c1=colorbar;
nexttile
h2 = imagesc(maskedOrientation); colormap turbo; axis image; axis off; title('Orientation')
set(h2, 'AlphaData', ~isnan(maskedOrientation));
c2=colorbar;
c2.Label.String = 'Degrees';
%% Hough transform

tic
sHough = GWHT(maskedOrientation,maskedMagnitude);
toc
%% 
% 

figure();
imagesc(sHough.hough, 'XData',sHough.theta, 'YData',sHough.rho);
colormap hsv; colorbar
xlabel('\theta'); ylabel('\rho'); 
title('GWHT')
%%
figure('Position',[294.6000 338 753.4000 420]);
imagesc(sHough.hough, 'XData',sHough.theta, 'YData',sHough.rho);
colormap turbo;
xlabel('\theta'); ylabel('\rho'); 
title('GWHT')
xline(-180,'g--')
xline(+180,'g--')
c=colorbar;
c.Label.String = 'Accumulator value';
%% Extend Hough space
% Hough space is now inherently extended to [-210,210] degrees within <./GWHT.m 
% GWHT.m>.
% 
% Peak relocating and clustering will happen after peaks are found.
%% Finding peaks

integerHough = round(sHough.hough.*100); %Multiplication to improve accuracy of rounding
rhoSpace = sHough.rho;
numPeaks = 100;

if useISAR
    peakThresh = 0.05*max(integerHough(:));
    peaks = houghpeaks(integerHough, ...
    numPeaks, ...
    "Threshold", ceil(peakThresh),...
    "NHoodSize", [501,51]);
else
    peakThresh = 0.1*max(integerHough(:));
    peaks = houghpeaks(integerHough, ...
    numPeaks, ...
    "Threshold", ceil(peakThresh))
end
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

figure('Position',[294.6000 338 753.4000 420]);
imagesc(sHough.hough, 'XData',sHough.theta, 'YData',sHough.rho);
colormap turbo;
c=colorbar;
c.Label.String = 'Accumulator value';
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

minPts = 1;

if useISAR
    % Need to do normalisation, or you'd think mahalanobis would account
    % for this
    eps = 0.1; % Initial test
    clusterLabel = dbscan(peaks(:,1:2), eps, minPts, 'Distance','mahalanobis');
else
    eps = 90;
    clusterLabel = dbscan(peaks, eps, minPts);
end
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
figure('Position',[294.6000 338 753.4000 420]);
gscatter(sHough.theta(peaks(:,2)), sHough.rho(peaks(:,1)), clusterLabel)
axis ij
ylim([min(sHough.rho), max(sHough.rho)])
xlim([-210 210])
xlabel('\theta'); ylabel('\rho'); 
title('Clustered peaks')
c=colorbar;
c.Label.String = 'Clusters';
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
for i=1:height(peaks)
    if peaks(i,2)<stopIdx
        corrTheta = round(peaks(i,2) + indexDiff); % This has to be an integer
        corrPoint = [peaks(i,1), corrTheta];

        % Check if this exists
        idx = find(ismember(peaks(:,1:2), corrPoint, "rows"));
        if ~isempty(idx)
            % Correlate these points
            mergingPairs = [mergingPairs; clusterLabel(i), clusterLabel(idx)];
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

figure('Position',[294.6000 338 753.4000 420]);
imagesc(sHough.hough, 'XData',sHough.theta, 'YData',sHough.rho);
colormap turbo;
c=colorbar;
c.Label.String = 'Accumulator value';
xlabel('\theta'); ylabel('\rho'); 
title('GWHT (w/ single peaks)')

x = min(sHough.theta) + meanPeaks(:,2)*thetaRes;
y = min(sHough.rho) + meanPeaks(:,1)*rhoRes;
% x = sHough.theta(peaks(:,2)); y = sHough.rho(peaks(:,1));
hold on
plot(x,y,'s','color','white');
hold off
xline(-150, 'g--')
%% Reconstruct lines

if useISAR
    meanPeaks = peaks(:,1:2);
end
for j=1:length(meanPeaks)
    pk = round(meanPeaks(j,:));
    % Find a new mask for that peak. Within a range of the peak angle
    pkAngle = wrapTo180(sHough.theta(pk(2))); % Need to wrap this back into -180:180 range
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
colours = hsv(length(uniqueVals));

for i=1:length(uniqueVals) %for each unique value
    idx = and(any(thetaRhoVals == uniqueVals(i,1),2), ...
        any(thetaRhoVals == uniqueVals(i,2),2)); % Matching rows
    tLines.label(idx) = i;
end

figure('Position',[295.4000 464.2000 754.4000 250.4000]);
if useISAR
    imshow(image,[])
else
    imshow(Inoise)
end
hold on
%% 
% Want it to be colour-coded to the cluster

for k = 1:height(tLines)
   xy = [tLines.point1(k,:); tLines.point2(k,:)];
   colour = colours(tLines.label(k),:);
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color',colour);
end
title('Line detection overlay');

hold off
%% Infinite lines
% Plot the infinite lines on the binary significance map.

whiteBackground = ones(size(mask));
for j=1:length(meanPeaks)
    pk = round(meanPeaks(j,:));
    % Find a new mask for that peak. Within a range of the peak angle
    pkAngle = wrapTo180(sHough.theta(pk(2))); % Need to wrap this back into -180:180 range
    pkLine = houghlines(whiteBackground, sHough.theta,rhoSpace,pk,...
        "FillGap",50, "MinLength", 100);

    if j==1 % For the first one, initialise the new big struct
        sLines = pkLine;
    else
        sLines = [sLines, pkLine];
    end
end
tLines = struct2table(sLines);
thetaRhoVals = [sLines.theta;sLines.rho]';
uniqueVals = unique(thetaRhoVals,"rows");
colours = hsv(length(uniqueVals));

for i=1:length(uniqueVals) %for each unique value
    idx = and(any(thetaRhoVals == uniqueVals(i,1),2), ...
        any(thetaRhoVals == uniqueVals(i,2),2)); % Matching rows
    tLines.label(idx) = i;
end

figure();
imshow(mask)
hold on
%% 
% Want it to be colour-coded to the cluster

for k = 1:height(tLines)
   xy = [tLines.point1(k,:); tLines.point2(k,:)];
   colour = colours(tLines.label(k),:);
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color',colour);
end
title('Line detection overlay');

hold off