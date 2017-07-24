This repository contains the code for all five hands-on of the ENS MVA 2014 lecture Vision 3D Artificielle of Mr Renaud Marlet and Mr Pascal Monasse (http://imagine.enpc.fr/~monasse/Stereo/)

Each folder contains the code and report for said hands-on.

Note that there is some mistake in the TP1 explaining the color mismatch.

TP1 : Panorama reconstruction by finding the (affine) transformation that produced the change in view in the images. The full panorama is reconstructed by doing a pull technique

![Panorama reconstruction](http://www.alth.fr/gallery/TP_vision/TP1_Full_compressed.jpg)

TP2 : We use a RANSAC algorithm to estimate the fundamental matrix enabling to relate the two images.

![Epipolar geometry](http://www.alth.fr/gallery/TP_vision/TP2_2_compressed.jpg)

TP3 : Harris algorithm to detect corner points is implemented. We also evaluate the adaptive non maxima suppression algorithm

![Harris corner detection](http://www.alth.fr/gallery/TP_vision/TP3_compressed.jpg)

TP4 : We implement the computation of disparity maps using seed propagation

![Disparity maps by seed propagation](http://www.alth.fr/gallery/TP_vision/TP4_12_compressed.jpg)

TP5 : We compute the disparity maps using graph-cut method. Additionally, graph-cut is used for segmentation.

![Disparity maps by graph cuts](http://www.alth.fr/gallery/TP_vision/TP5_12_compressed.jpg)
