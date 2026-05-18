function sHoughGW = GWHT(maskedOrientation, maskedMagnitude, thetaRes, sigma)
%GWHT  Gradient Weighted Hough Transform
%
% DOCSTRING
%   Computes an image-gradient-weighted Hough transform, using pixel
%   gradient magnitude and orientation to weight votes in the accumulator
%   array. Orientation and magnitude images should be masked to significant
%   pixels beforehand to reduce computation time and suppress noise.
% 
%   Returns a struct containing the accumulator array along with the
%   corresponding theta and rho axes
%
%   Inputs:
%       maskedOrientation — 2‑D array of pixel orientation values (degrees)
%       maskedMagnitude   — 2‑D array of pixel gradient magnitudes
%       thetaRes          — Angular resolution in degrees (default: 0.2)
%       sigma             — Standard deviation of the Gaussian angular kernel 
%                           (default: 20)
%
%   Output:
%       sHoughAngular — Struct containing:
%           .hough    — Angular Hough accumulator
%           .theta    — Theta axis values (deg)
%           .rho      — Rho axis values (pixels)

arguments
    maskedOrientation (:,:) double
    maskedMagnitude   (:,:) double

    thetaRes (1,1) double {mustBePositive} = 0.2
        % Angular resolution in degrees

    sigma (1,1) double {mustBePositive} = 10
        % Standard deviation of Gaussian angular kernel
end

%% Input errors and warnings

% Ensure images are the same size
if ~isequal(size(maskedOrientation),size(maskedMagnitude))
    error('Input image arrays must be the same size. Currently %s and %s.', mat2str(size(maskedMagnitude)), mat2str(size(maskedOrientation)))
end

% Check if images have been masked
if ~any(isnan(maskedOrientation(:))) || ~any(isnan(maskedMagnitude(:)))
    warning('One or both input images have not been masked. Performance will be signficantly reduced.')
    userResponse = input("Input 'Y' to continue, or anything else to exit: ",'s');
    if ~strcmpi(userResponse, 'y') % Capital and lowercase both work
        return;
    end
end
%% Helper variable initialisation

thetaVals = -210:thetaRes:210;
[x, y] = size(maskedOrientation);
rhoMaximum = norm([x y]);
rhoSpace = (-rhoMaximum:1:rhoMaximum);

% Number of theta and rho bins
num_thetas = numel(thetaVals);
num_rhos = numel(rhoSpace);

% Accumulator array
houghSpace = zeros(num_rhos, num_thetas);

% Define gaussian function
gaussFunc = circshift(normpdf(thetaVals,90,sigma),-90/thetaRes); % Sigma should be tuned
%% Gradient weighted Hough transform calculation

for xi = 1:x
    for yj = 1:y
        % For each non-NaN pixel
        if ~isnan(maskedOrientation(xi, yj))
            pixelOrientation = round(maskedOrientation(xi,yj)/thetaRes)*thetaRes;

            pixelMagnitude = maskedMagnitude(xi,yj);

            pixelGauss = circshift(gaussFunc,round(pixelOrientation/thetaRes)); % Shift the gauss function to be over the orientation of the pixel
            for thetaIdx = 1:num_thetas
                % For each theta, calculate rho at that (x,y) position
                theta = thetaVals(thetaIdx);

                % Parametric representation of a line
                rho = yj * cosd(theta) + xi * sind(theta);
                
                % Find corresponding index in rhoSpace
                rhoIdx = round(rho + rhoMaximum + 1);
                
                % Calculate vote strength and add to accumulator
                pixelVote = pixelGauss(thetaIdx)*pixelMagnitude;
                houghSpace(rhoIdx, thetaIdx) = ...
                    houghSpace(rhoIdx, thetaIdx) + pixelVote;
            end
        end
    end
end

%% Outputs

sHoughGW.hough = houghSpace;
sHoughGW.theta = thetaVals;
sHoughGW.rho = rhoSpace;
end