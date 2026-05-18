function [Gxa, Gya] = gradientByRatio_v2(image, alpha, gridSize)
%GRADIENTBYRATIO_V2  Compute gradient by ratio of an image.
%
%   [Gxa, Gya] = GRADIENTBYRATIO_V2(image, alpha, gridSize) computes 
%   horizontal and vertical gradient estimates by convolving the image with 
%   paired anisotropic exponential filters and taking the log‑ratio of their 
%   responses. Filters are generated from an L1-based exponential kernel and 
%   cached using persistent variables to avoid recomputation when alpha and 
%   gridSize remain unchanged.
%
%   INPUTS:
%       image      - Numeric 2D image array.
%       alpha      - Positive scalar controlling the exponential decay.
%       gridSize   - Integer specifying half-size of the filter window 
%                    (default = 20). Full filter size is (2*gridSize+1).
%
%   OUTPUTS:
%       Gxa        - Log‑ratio gradient estimate in the horizontal direction.
%       Gya        - Log‑ratio gradient estimate in the vertical direction.

arguments
    image
    alpha (1,1) double {mustBePositive}
    gridSize (1,1) double {mustBeInteger} = 20
end

% Persistent cached filters
persistent lastAlpha lastGrid fExpLeft fExpRight fExpTop fExpBottom

% Rebuild filters if alpha or gridsize has changed
if isempty(lastAlpha) || isempty(lastGrid) || lastAlpha ~= alpha || lastGrid ~= gridSize

    % Store new settings
    lastAlpha = alpha;
    lastGrid = gridSize;

    [X,Y] = meshgrid(-gridSize:gridSize, -gridSize:gridSize);
    exponentialFilter = exp(-(abs(X)+abs(Y))./alpha);
    
    maskTop = zeros(2*gridSize+1);
    maskTop(1:gridSize, :) = 1;
    fExpTop = maskTop .* exponentialFilter;
    
    % Define other filters
    fExpBottom = flipud(fExpTop);
    fExpLeft = fExpTop';
    fExpRight = fliplr(fExpLeft);
    
    % Normalisation 
    fExpTop = fExpTop ./ sum(fExpTop, "all");
    fExpBottom = fExpBottom ./ sum(fExpBottom, "all");
    fExpLeft = fExpLeft ./ sum(fExpLeft, "all");
    fExpRight = fExpRight ./ sum(fExpRight, "all");
end
R1 = imfilter(image, fExpLeft, 'replicate')./imfilter(image, fExpRight, 'replicate');
R3 = imfilter(image, fExpTop, 'replicate')./imfilter(image, fExpBottom, 'replicate');

Gxa = log(R1);
Gya = log(R3);
end