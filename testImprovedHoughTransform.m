% test images
addpath('TestData\')
 % = imread("TestData\testImage.png");

% Test data files:
% - Original image (.png)
% - gradientByRatio_INPUT (.mat)
% - MagnitudeOrientationImages (.mat)
% Check that 

%% Check for required toolboxes
[files, requiredAddons] = matlab.codetools.requiredFilesAndProducts('ImprovedHoughTransform.m');

tInstalledAddons = matlab.addons.installedAddons;

for i=1:length(requiredAddons)
    if ~strcmp(requiredAddons(i).Name, 'MATLAB') % 'MATLAB' doesn't appear in the addons list.
        assert(any(strcmp(tInstalledAddons.Name,requiredAddons(i).Name)), ...
            ['Missing addon: ' char(requiredAddons(i).Name)])
    end
end

%% Testing gradientByRatio_v2
% Check that applying gradientByRatio_v2 to gradientByRatio_INPUT yields
% magnitude and orientation images as expected
inputImage = load('TestData\Exemplary lines\gradientByRatio_INPUT.mat').Inoise;
load('TestData\Exemplary lines\calculatedGradients.mat')

%% GR: Test valid inputs
alpha = 20;
gridSize = 10;

[Gx, Gy] = gradientByRatio_v2(inputImage, alpha, gridSize);
assert(~isempty(Gx))
assert(~isempty(Gy))

%% GR: Test gridSize default behaviour
alpha = 10;

[outputGxa, outputGya] = gradientByRatio_v2(inputImage, alpha);
assert(~isempty(outputGxa))
assert(~isempty(outputGya))

%% GR: Test expected output (using known input/output data)
[outputGxa, outputGya] = gradientByRatio_v2(inputImage, 20);
assert(isequal(Gxa,outputGxa))
assert(isequal(Gya,outputGya))

%% Testing GWHT

% Check that the output is as expected

% Check that passing bad data will yield the errors and warnings as
% expected.