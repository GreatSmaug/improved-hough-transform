# improved-hough-transform
Improved Hough Transform workflow, including gradient-weighted Hough transform (GWHT) function and robust edge detection methodology for SAR/ISAR imagery.

## _ImprovedHoughTransform.m_
Main script file. Live script version available as .mlx.

## _StandardHoughTransform.m_
Implementation of the standard Hough transform, included as a contrast to the improved methodology.

## Supporting files
### TestData directory
Contains the data needed to run the testfile _testImprovedHoughTransform.m_ 
```
TestData/
└── Exemplary Lines/
    ├── calculatedGradients.mat
    ├── gradientByRatio_INPUT.mat
    └── originalImage.png
```

### testImage.png
Default image to be used by _ImprovedHoughTransform.m_

### ISARexample.mat
Example simulated ISAR image, which can be used in _ImprovedHoughTransform.m_ or _StandardHoughTransform.m_ by setting _useISAR=true_.

### _testImprovedHoughTransform.m_
Testfile containing tests for gradientByRatio_v2 function. Can be run from MATLAB console by:
```
runtests
```
Tests for GWHT will be added in a future update.

## Functions
### gradientByRatio_v2.m
Implementation of gradient-by-ratio image gradient calculation, recommended for SAR/ISAR imagery due to speckle noise [1,2].

### GWHT.m
The **gradient-weighted Hough transform**, an adaptation of the standard Hough transform using calculated image gradients to weight accumulator voting. Aims to reduce noise and artefacting in Hough space, allowing subtler peaks to be detected more easily. Methodology presented in my own paper [3].



## References
[1] R. Fjortoft, A. Lopes, P. Marthon, and E. Cubero-Castan, ‘An optimal multiedge detector for SAR image segmentation’, IEEE Transactions on Geoscience and Remote Sensing, vol. 36, no. 3, pp. 793–802, May 1998, doi: 10.1109/36.673672.

[2] F. Dellinger, J. Delon, Y. Gousseau, J. Michel, and F. Tupin, ‘SAR-SIFT: A SIFT-Like Algorithm for SAR Images’, IEEE Transactions on Geoscience and Remote Sensing, vol. 53, no. 1, pp. 453–466, Jan. 2015, doi: 10.1109/TGRS.2014.2323552.

[3] M. Coe, G. Jones, L.-N. Alconcel, and M. Gashinova, ‘Persistent feature reconstruction of resident space objects (RSOs) within inverse synthetic aperture radar (ISAR) images’, Dec. 17, 2025, arXiv: arXiv:2512.15618. doi: 10.48550/arXiv.2512.15618.

