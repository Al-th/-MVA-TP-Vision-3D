// Imagine++ project
// Project:  HarrisANMS
// Author:   Alexandre THIS
// Date:     2014/10/15

#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <Imagine/Images.h>
#include <limits>
#include <vector>
#include <algorithm>
#include <cstdlib>

#define HARRIS_RESPONSE 0
#define HARRIS_STRENGTH 1

#define STANDARD_NMS 0
#define ADAPTATIVE_NMS 1

#define VERBOSE 1

using namespace Imagine;
using namespace std;

typedef struct responsePoint{
  int x;
  int y;
  float r;
  float radius;
}t_responsePoint;


Vector<float> gaussianKernel(float sigma){
  int n = round(3*sigma)-1;
  if(n<1){
    cout << "Exiting software, bad gaussian kernel allocation, kernel size is too small"<< endl;
    exit(0);
  }
  int sum = 0;
  Vector<float> g(2*n+1);
  for(int i = -n; i < n+1; i++){
    g[i+n] = exp(-(i*i)/(2*sigma*sigma));
    sum += g[i+n];
  }
  g/=sum;
  return g;
}

Vector<float> padWithZeros(const Vector<float> v, int zeroPaddings){
  int s = v.size();
  s += 2*zeroPaddings;
  Vector<float> newV(s);
  for(int i = 0; i < s; i++){
    newV[i] = 0;
  }
  
  for(int i = 2; i < s-2;i++){
    newV[i] = v[i-2];
  }
  return newV;
  
}

void displayMatrix(const Matrix<float> Mat){
  Matrix<float> Mat2(Mat.width(),Mat.height());
  Image<Color,2> IRes(Mat.width(),Mat.height());
  float min = numeric_limits<float>::max();
  float max = numeric_limits<float>::min();
  
  for(int i = 0; i < Mat.width();i++){
    for(int j = 0; j < Mat.height(); j++){
      if(Mat(i,j) < min){
	min = Mat(i,j);
      }
      if(Mat(i,j) > max){
	max = Mat(i,j);
      }
      Mat2(i,j) = Mat(i,j);
    }
  }
  
   for(int i = 0; i < Mat.width();i++){
    for(int j = 0; j < Mat.height(); j++){
      Mat2(i,j) -= min;
    }
   }
   Mat2 /= (max-min);
   Mat2 *= 255;
  
  
  for(int i = 0; i < Mat.width(); i++){
    for(int j = 0; j < Mat.height(); j++){
      IRes(i,j) =  Color(Mat2(i,j),Mat2(i,j),Mat2(i,j));
    }
  }
  Window w = openWindow(IRes.width(), IRes.height());
  setActiveWindow(w);
  display(IRes,0,0);
}

Vector<float> convolve1D1D(const Vector<float> v1, const Vector<float> kernel){
  Vector<float> convolved(v1.size());
  for(int i = 0; i < convolved.size();i++){
    convolved[i] = 0;
  }
  if(kernel.size()%2!=1){
    cout << "Problem with kernel, size is not odd" << endl;
  }
  else{
    float sumKernel = 0;
    for(int i = 0; i < kernel.size(); i++){
      sumKernel += kernel[i];
    }
    int boundaryOffset = floor(kernel.size()/2);
    for(int i = 0; i < boundaryOffset;i++){
      convolved[i] = 0;
    }
    for(int i = boundaryOffset; i < v1.size()-boundaryOffset;i++){
      convolved[i] = 0;
      for(int k = 0; k < kernel.size();k++){
	convolved[i] += kernel[k]*v1[i+k-boundaryOffset];
      }
      if(sumKernel != 0){
	convolved[i] / sumKernel;
      }
    }
  }
  return convolved;
}

Matrix<float> convolve2D1D(const Matrix<float> I, const Vector<float> kernel, int dim){
  Matrix<float> convolved(I.width(),I.height());
  for(int i = 0; i < I.width(); i++){
    for(int j = 0; j < I.height(); j++){
      convolved(i,j) = 0;
    }
  }
  float minVal = 2048;
  float maxVal = 0;
  if(kernel.size()%2!=1){
    cout << "Problem with kernel, size is not odd" << endl;
    return convolved;
  }
  else{
    float sumKernel = 0;
    for(int i = 0; i < kernel.size(); i++){
      sumKernel += kernel[i];
    }
    int boundaryOffset = floor(kernel.size()/2);
    for(int i = ((dim==0)*boundaryOffset); i < I.width()-((dim==0)*boundaryOffset); i++){
      for(int j = ((dim==1)*boundaryOffset); j < I.height()-((dim==1)*boundaryOffset); j++){
	convolved(i,j) = 0;
	for(int k = 0; k < kernel.size(); k++){
	  float kernelValue = kernel[k];
	  convolved(i,j) += kernelValue*I(i+((dim==0)*(k-boundaryOffset)),j+((dim==1)*(k-boundaryOffset)));
	}
	if(sumKernel != 0){
	  convolved(i,j) /= sumKernel;
	}
      }
    }
  }
  return convolved;
}

bool responseComparator(const t_responsePoint& p1, const t_responsePoint& p2){
  return (p1.r > p2.r); 
}

bool radiusComparator(const t_responsePoint& p1, const t_responsePoint& p2){
  return(p1.radius > p2.radius);
}



void standardNonMaximaSuppresion(const Image<Color,2> I1, const Matrix<float> HarrisResponse, int nbPoints, char* fileName){
  vector<t_responsePoint> responsePoint;
  for(int i = 1; i < I1.width()-1; i++){
    for(int j = 1; j < I1.height()-1;j++){
      t_responsePoint p;
      p.x = i;
      p.y = j;
      p.r = HarrisResponse(i,j);
      bool isMax = true;
      for(int k = -1; k < 2; k++){
	for(int l = -1; l < 2; l++){
	  if(p.r < HarrisResponse(i+k,j+l)){
	    isMax = false;
	  }
	}
      }
      if(isMax){
	responsePoint.push_back(p);
      }
    }
  }

  sort(responsePoint.begin(), responsePoint.end(), responseComparator);
  
  Window w = openWindow(I1.width(), I1.height(), fileName);
  setActiveWindow(w);
  display(I1,0,0);
  
  for(int i = 0; i < nbPoints ;i++){
    fillCircle(responsePoint.at(i).x, responsePoint.at(i).y,3,RED);
  }
}

void adaptativeNonMaximaSuppression(const Image<Color,2> I1, const Matrix<float> HarrisResponse, int nbPoints, char* fileName){
  vector<t_responsePoint> sortedResponsePoint;
  for(int i = 0; i < I1.width(); i++){
    for(int j = 0; j < I1.height(); j++){
      t_responsePoint p;
      p.x = i;
      p.y = j;
      p.r = HarrisResponse(i,j);
      sortedResponsePoint.push_back(p);
    }
  }
  
  
  sort(sortedResponsePoint.begin(), sortedResponsePoint.end(), responseComparator);
  
  t_responsePoint p = sortedResponsePoint.at(0);
  p.radius = 1000000;

  vector<t_responsePoint> radiusPoints;
  radiusPoints.push_back(p);


  for(int j = 0; j < sortedResponsePoint.size();j++){
    t_responsePoint p = sortedResponsePoint[j];
    float radius = 10000000; // big value
    for(int i = 0; i < radiusPoints.size(); i++){
      t_responsePoint q = radiusPoints[i];
      if(p.r < 0.9*q.r){
	float currRadius = (p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y);
	if(currRadius < radius){
	  radius = currRadius;
	}
      }
    }
    if(radius == 10000000){
      radius = 0;
    }
    p.radius = radius;
    radiusPoints.push_back(p);
   
  }


  sort(radiusPoints.begin(), radiusPoints.end(),radiusComparator);
  
  Window w = openWindow(I1.width(), I1.height(), fileName);
  setActiveWindow(w);
  display(I1,0,0);
  
  for(int i = 0; i < nbPoints ;i++){
    fillCircle(radiusPoints.at(i).x, radiusPoints.at(i).y,3,RED);
  }
}

void computeCorners(const Image<Color,2> I1, int typeDetector, int typeNMS, int nbPoints, float sigma_d, float sigma_i){
  float lambdaResponse = 0.04 ;
   
  char* outputFileName = (char*)malloc(256*sizeof(char));
  sprintf(outputFileName, "%d-%d-%d-%2.2f-%2.2f", typeDetector,typeNMS,nbPoints,sigma_d,sigma_i);
    


  Matrix<float> I_Matrix(I1.width(),I1.height());
  for(int i = 0; i < I1.width(); i++){
    for(int j = 0; j < I1.height(); j++){
      I_Matrix(i,j) = I1(i,j);
    }
  }
  
  Vector<float> derivationOperator(3);
  derivationOperator[0] = -1;
  derivationOperator[1] = 0;
  derivationOperator[2] = 1;
  
  Vector<float> gNonPadded = gaussianKernel(sigma_d);
  Vector<float> gaussian = padWithZeros(gNonPadded,2);

  Vector<float> dx = convolve1D1D(gaussian,derivationOperator);
  
  Matrix<float> Ix(I1.width(),I1.height()),Iy(I1.width(),I1.height());
  Ix = convolve2D1D(I_Matrix,dx,0);
  Iy = convolve2D1D(I_Matrix,dx,1);

  Vector<float> g2 = gaussianKernel(sigma_i);
  
  Ix = convolve2D1D(Ix,g2,0);
  Iy = convolve2D1D(Iy,g2,1);
  
  

  Matrix<float> Ixx(I1.width(),I1.height());
  Matrix<float> Ixy(I1.width(),I1.height());
  Matrix<float> Iyy(I1.width(),I1.height());  

  for(int i = 0; i < I1.width();i++){
    for(int j = 0; j < I1.height(); j++){
      Ixx(i,j) = Ix(i,j)*Ix(i,j);
      Iyy(i,j) = Iy(i,j)*Iy(i,j);
      Ixy(i,j) = Ix(i,j)*Iy(i,j);
      
      
    }
  }
  
  Ixx = convolve2D1D(Ixx,g2,0);
  Ixx = convolve2D1D(Ixx,g2,1);
  Iyy = convolve2D1D(Iyy,g2,0);
  Iyy = convolve2D1D(Iyy,g2,1);
  Ixy = convolve2D1D(Ixy,g2,0);
  Ixy = convolve2D1D(Ixy,g2,1);

  Matrix<float> HarrisResponse(I1.width(),I1.height());
  


  for(int i = 0 ; i < I1.width(); i++){
    for(int j = 0; j < I1.height(); j++){
      FMatrix<float,2,2> A(0.0);
      A(0,0) = Ixx(i,j);
      A(1,1) = Iyy(i,j);
      A(1,0) = Ixy(i,j);
      A(0,1) = Ixy(i,j);
      float detA = A(0,0)*A(1,1)-A(1,0)*A(0,1);
      float trA = A(0,0)+A(1,1);
      
      if(typeDetector == HARRIS_STRENGTH){
	if(trA != 0){
	HarrisResponse(i,j) = detA / trA;
	}
	else{
	  HarrisResponse(i,j) = 0;
	}
      }
      else{
	HarrisResponse(i,j) = detA - lambdaResponse*(trA*trA);
      }
    }
  }

  if(VERBOSE){
    displayMatrix(Ix);
    displayMatrix(Iy);
    displayMatrix(Ixx);
    displayMatrix(Ixy);
    displayMatrix(Iyy);
    displayMatrix(HarrisResponse);
  }
  if(typeNMS == STANDARD_NMS){
    standardNonMaximaSuppresion(I1,HarrisResponse,nbPoints, outputFileName);
  }
  else{
    adaptativeNonMaximaSuppression(I1,HarrisResponse,nbPoints, outputFileName);
  }
  
}


int main(int argc, char* argv[])
{
    const char* s1 = argc>1? argv[1]: srcPath("corners.jpg");

    // Load and display images
    Image<Color,2> I1;
    if( ! load(I1, s1) ) {
        cerr<< "Unable to load image" << endl;
        return 1;
    }
   
    computeCorners(I1, HARRIS_RESPONSE,STANDARD_NMS, 65,0.5,1);

    endGraphics();
    return 0;
    
}
