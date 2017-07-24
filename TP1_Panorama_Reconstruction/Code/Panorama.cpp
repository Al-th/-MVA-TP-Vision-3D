// Imagine++ project
// Project:  Panorama
// Author:   Pascal Monasse
// Date:     2012/09/25

#include <Imagine/Graphics.h>
#include <Imagine/Images.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <sstream>

using namespace Imagine;
using namespace std;

// Record clicks in two images, until right button click
void getClicks(Window w1, Window w2,
               vector<IntPoint2>& pts1, vector<IntPoint2>& pts2) {

  //Variables used to recover informations about clicks in the image
  IntPoint2 p;
  int button;

  //Show the first window to the user and wait for him to perform clicks action. If click is not a Right click, add an interest point and show it to the user thanks to a circle and a number.
  setActiveWindow(w1);
  showWindow(w1);
  cout << "Selectionnez les points d'interets de la premiere image avec le clic gauche. Appuyez sur le clic droit pour terminer." << endl;
  cout << endl;
  while((button=getMouse(p))!=3){
    pts1.push_back(p);
    fillCircle(p,5,RED);
    char* nbElts = new char[5];
    sprintf(nbElts,"%d",pts1.size());
    drawString(p,nbElts,BLUE,17,0,0,true,false);
  }
  
  //Perform the same task on the second window.
  setActiveWindow(w2);
  showWindow(w2); 
  cout << "Selectionnez les points d'interets correspondants de la seconde image avec le clic gauche. Appuyez sur le clic droit pour terminer." << endl;
  cout << endl;
  while((button=getMouse(p))!=3){
    pts2.push_back(p);
    fillCircle(p,5,RED);
    char* nbElts = new char[5];
    sprintf(nbElts,"%d",pts2.size());
    drawString(p,nbElts,BLUE,17,0,0,true,false);
  }

  
}

// Return homography compatible with point matches
Matrix<float> getHomography(const vector<IntPoint2>& pts1,
                            const vector<IntPoint2>& pts2) {
    size_t n = min(pts1.size(), pts2.size());
    if(n<4) {
        cout << "Not enough correspondences: " << n << endl;
        return Matrix<float>::Identity(3);
    }
    Matrix<float> A(2*n,8);
    Vector<float> B(2*n);
    
    //Filling the matrix A and the vector B to prepare the linear system to be solved (Ah=B)
    for(int i = 0; i < n; i++){
      float x,y,xp,yp;
      
      x = (float)(pts2[i].x());
      y = (float)(pts2[i].y());
      xp = (float)(pts1[i].x());
      yp = (float)(pts1[i].y());
      
      
      A(2*i,0) = x;
      A(2*i,1) = y;
      A(2*i,2) = 1;
      A(2*i,3) = 0;
      A(2*i,4) = 0;
      A(2*i,5) = 0;
      A(2*i,6) = -(xp*x);
      A(2*i,7) = -(xp*y);
      B[2*i] = xp;

      A(2*i+1,0) = 0;
      A(2*i+1,1) = 0;
      A(2*i+1,2) = 0;
      A(2*i+1,3) = x;
      A(2*i+1,4) = y;
      A(2*i+1,5) = 1;
      A(2*i+1,6) = -(yp*x);
      A(2*i+1,7) = -(yp*y);
      B[2*i+1] = yp;
    }
   
    B = linSolve(A, B);
    Matrix<float> H(3, 3);

    H(0,0)=B[0]; H(0,1)=B[1]; H(0,2)=B[2];
    H(1,0)=B[3]; H(1,1)=B[4]; H(1,2)=B[5];
    H(2,0)=B[6]; H(2,1)=B[7]; H(2,2)=1;

    // Sanity check
    cout << "Sanity Check" << endl;
    for(size_t i=0; i<n; i++) {
        float v1[]={(float)pts1[i].x(), (float)pts1[i].y(), 1.0f};
        float v2[]={(float)pts2[i].x(), (float)pts2[i].y(), 1.0f};
        Vector<float> x1(v1,3);
        Vector<float> x2(v2,3);
        x1 = H*x1;
        cout << x1[1]*x2[2]-x1[2]*x2[1] << ' '
             << x1[2]*x2[0]-x1[0]*x2[2] << ' '
             << x1[0]*x2[1]-x1[1]*x2[0] << endl;
    }
    cout << endl;
    return H;
}

// Grow rectangle of corners (x0,y0) and (x1,y1) to include (x,y)
void growTo(float& x0, float& y0, float& x1, float& y1, float x, float y) {
    if(x<x0) x0=x;
    if(x>x1) x1=x;
    if(y<y0) y0=y;
    if(y>y1) y1=y;    
}

// Panorama construction
void panorama(const Image<Color,2>& I1, const Image<Color,2>& I2,
              Matrix<float> H) {
    Vector<float> v(3);
    float x0=0, y0=0, x1=I2.width(), y1=I2.height();

    v[0]=0; v[1]=0; v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=I1.width(); v[1]=0; v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=I1.width(); v[1]=I1.height(); v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=0; v[1]=I1.height(); v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    cout << "x0 x1 y0 y1=" << x0 << ' ' << x1 << ' ' << y0 << ' ' << y1<<endl;

    Image<Color> I(int(x1-x0), int(y1-y0));
    setActiveWindow( openWindow(I.width(), I.height()) );
    I.fill(WHITE);


    //Panorama reconstruction with "Pull" method

    //Inverse the homography transformation
    Matrix<float> inversedHomography(3,3);
    inversedHomography = inverse(H);

    //For each pixels in the final images
    for(int i = 0; i < I.width();i++){
      for(int j = 0; j < I.height();j++){
	//We compute the position of the point (i,j) before the homography thanks to the inverse homography 
	//We also normalize to get homogeneous coordinates and round its values to the nearest integer (It is indeed not the best way to interpolate the pixel color)
	Vector<float> p(3);
	p[0] = (float)i; 
	p[1] = (float)j; 
	p[2] = 1.0f;
	p = inversedHomography*p;
	p /= p[2];
	for(int k = 0; k < 3;k++){
	  p[k] = round(p[k]);
	}

	//We fill the image I with I1 and I2. If there is a superposition between the two images I1 and I2, we average the color of the pixel at the position (i,j)
	bool pixInI1 = false;
	bool pixInI2 = false;
	if(p[0] >= 0 && p[0] < I2.width() && p[1] >= 0 && p[1] < I2.height()){
	  I(i,j) = I2(p[0],p[1]);
	    pixInI2 = true;
	  
	}
	if(i >= 0 && i < I1.width() && j >= 0 && j < I1.height()){
	  I(i,j) += I1(i,j);
	    pixInI1 = true;
	}
	
	if(pixInI1 && pixInI2){
	  I(i,j) /= 2;
	}
      }
    }
    
    
    display(I,0,0);
}

// Main function
int main() {
    const char* s1 = srcPath("image0006.jpg");
    const char* s2 = srcPath("image0007.jpg");

    // Load and display images
    Image<Color> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load the images" << endl;
        return 1;
    }
    Window w1 = openWindow(I1.width(), I1.height(), s1);
    display(I1,0,0);
    Window w2 = openWindow(I2.width(), I2.height(), s2);
    setActiveWindow(w2);
    display(I2,0,0);

    // Get user's clicks in images
    vector<IntPoint2> pts1, pts2;
    getClicks(w1, w2, pts1, pts2);

    vector<IntPoint2>::const_iterator it;
    cout << "pts1="<<endl;
    for(it=pts1.begin(); it != pts1.end(); it++)
        cout << *it << endl;
    cout << "pts2="<<endl;
    for(it=pts2.begin(); it != pts2.end(); it++)
        cout << *it << endl;

    // Compute homography
    Matrix<float> H = getHomography(pts1, pts2);
    cout << "H=" << H/H(2,2);

    // Apply homography
    panorama(I1, I2, H);

    endGraphics();
    return 0;
}
