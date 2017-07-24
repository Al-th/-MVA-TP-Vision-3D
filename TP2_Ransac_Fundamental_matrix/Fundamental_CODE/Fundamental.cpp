// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse
// Date:     2013/10/08

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace Imagine;
using namespace std;

struct Match {
    float x1, y1, x2, y2;
};

// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color,2> I1, Image<Color,2> I2,
              std::vector<Match>& matches) {
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0,0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(),0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0*100.0;
    for(size_t i=0; i < feats1.size(); i++) {
        SIFTDetector::Feature f1=feats1[i];
        for(size_t j=0; j < feats2.size(); j++) {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if(d < MAX_DISTANCE) {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}


/*Author : Alexandre This
 *FMatrix<float,9,9> getRandomAMatrix(const std::vector<Match> & matches)
 *
 *Select at random 8 matches and create the A matrix to be used in 8-point algorithm
 */
FMatrix<float,9,9> getRandomAMatrix(const std::vector<Match> & matches){    
  
  FMatrix<float,3,3> normalizationMatrix(0.0);
  normalizationMatrix(0,0) = 0.001;
  normalizationMatrix(1,1) = 0.001;
  normalizationMatrix(2,2) = 1;
  
  FMatrix<float,9,9> A(0.0);
  for(int i = 0; i < 8; i++){
    FVector<float,3> pl(0),pr(0);
    

    //Select matches at random;
    int randIndex = rand()%matches.size();
    
    pl[0] = matches.at(randIndex).x1;
    pl[1] = matches.at(randIndex).y1;
    pl[2] = 1;
    
    pr[0] = matches.at(randIndex).x2;
    pr[1] = matches.at(randIndex).y2;
    pr[2] = 1;
    
    //normalize points
    pl = normalizationMatrix*pl;
    pr = normalizationMatrix*pr;
    
    //fill the matrix
    A(i,0) = pl[0]*pr[0] ;
    A(i,1) = pl[0]*pr[1] ;
    A(i,2) = pl[0] ;
    A(i,3) = pl[1]*pr[0] ;
    A(i,4) = pl[1]*pr[1] ;
    A(i,5) = pl[1] ;
    A(i,6) = pr[0] ;
    A(i,7) = pr[1] ;
    A(i,8) = 1 ;
  }
  
  //Add a dumb line at the end
  for(int i = 0; i < 9;i++){
    A(8,i) = 0;
  }
  
  return A;
}


/*Author : Alexandre This
 *FMatrix<float,3,3> compute8PointsFundamentalEstimation(FMatrix<float,9,9> & A)
 *
 *Use the A matrix (9*9 matrix with a dumb line) to compute an estimation of the Fundamental matrix using the 8-point algorithm.
 */
FMatrix<float,3,3> compute8PointsFundamentalEstimation(FMatrix<float,9,9> & A){
  FMatrix<float,3,3> reformedF(0.0);
  
  FMatrix<float,9,9> Ua(0.0),Vta(0.0);
  FVector<float,9> Sa(0.0);
  
  FMatrix<float,3,3> badRankedFMatrix(0.0);
  FMatrix<float,3,3> Ubf(0.0),Vtbf(0.0);
  FVector<float,3> Sbf(0.0);
  
  FMatrix<float,3,3> normalizationMatrix(0.0);
  normalizationMatrix(0,0) = 0.001;
  normalizationMatrix(1,1) = 0.001;
  normalizationMatrix(2,2) = 1;
  
  //Solve the minimisation problem
  A = transpose(A)*A;
  
  svd(A,Ua,Sa,Vta);
  
  for(int i = 0; i < 9; i++){
    badRankedFMatrix(i/3,i%3) = Vta(8,i);
  }
  
  //Force the rank to 2
  svd(badRankedFMatrix,Ubf,Sbf,Vtbf);
  Sbf[2] = 0;
  
  reformedF = Ubf*(Diagonal(Sbf)*Vtbf);
  
  //Unnormalize matrix
  reformedF = transpose(normalizationMatrix)*(reformedF*normalizationMatrix);
  
  //Ensure that the norm is 1
  reformedF = reformedF/norm(reformedF);

  return reformedF;
  
}


/*Author : Alexandre This
 *vector<Match> computeInliersForMatrixEstimate(const FMatrix<float,3,3> FEstimate, const std::vector<Match> & matches, float thresholdDistanceRansac)
 *
 *Use the estimated FEstimate matrix to compute inliers. Inliers are defined as the set of matches that fit best the current transformation. 
 *Fit best is defined as : the orthogonal distance between a point and the epipolar line is less than a thresold thresholdDistanceRansac
 */
vector<Match> computeInliersForMatrixEstimate(const FMatrix<float,3,3> FEstimate, const std::vector<Match> & matches, float thresholdDistanceRansac){
  vector<Match> inliers;
  
  for(int i = 0; i < matches.size(); i++){
    FVector<float,3> pl(0),pr(0);
    
    pl[0] = matches.at(i).x1;
    pl[1] = matches.at(i).y1;
    pl[2] = 1;
    
    pr[0] = matches.at(i).x2;
    pr[1] = matches.at(i).y2;
    pr[2] = 1;
    
    
    //Compute the epipolar line
    FVector<float,3> line = FEstimate*pr;
    
    //Normalize epipolar line to ensure avoiding numerical problems
    line = line/(sqrt(line[0]*line[0]+line[1]*line[1]));

    //compute orthogonal distance to the epipolar line
    float distance = pl*line;
    distance *= distance;
    
    if(distance < thresholdDistanceRansac){
      inliers.push_back(matches.at(i));
    }
  }
  
  return inliers;
}


/*Author : Alexandre This
 *vector<Match> computeRansacInliers(const std::vector<Match> & matches)
 *
 *Compute the matches fitting the best the transformation using the 8-point algorithm to estimate the F matrix, and ransac to find the best estimation of F
 *
 */
vector<Match> computeRansacInliers(const std::vector<Match> & matches){
  int nbInterations = 500;
  float thresholdDistanceRansac = 1; // 1 pixel 
  
  vector<Match> finalInliers;
  
  for(int iter = 0; iter < nbInterations; iter++){
    vector<Match> currentInliers;
    
    FMatrix<float,9,9> A(0.0);
    FMatrix<float,3,3> fundamentalEstimate(0.0);
    
    //Get a random A Matrix 
    A = getRandomAMatrix(matches);
    
    //Compute the fundamental matrix estimate 
    fundamentalEstimate = compute8PointsFundamentalEstimation(A);
    
    //Get the inliers for this particular fundamental estimate
    currentInliers = computeInliersForMatrixEstimate(fundamentalEstimate, matches, thresholdDistanceRansac);
    
    //Check if this estimate yields better results
    int currentInliersSize = currentInliers.size();
    int finalInliersSize = finalInliers.size();
    if(currentInliersSize > finalInliersSize){
      cout << "Ransac new size : " << currentInliersSize << endl;
      finalInliers.clear();
      for(int i = 0; i < currentInliersSize; i++){
	finalInliers.push_back(currentInliers.back());
	currentInliers.pop_back();
      }
    }
  }
  
  return finalInliers;
  
}


/*Author : Alexandre This
 *FMatrix<float,3,3> refineFundamentalMatrix(std::vector<Match> & bestRansacMatches)
 *
 *Refine the fundamental matrix using a vector of matches that best fit the transformation
 *
 *Compute the least square solution to the following minimization problem :
 *Min ||Af||2 subject to ||f|| = 1
 *
 *The solution is given by the eigenvector corresponding to the smallest eigenvalue of A
 */
FMatrix<float,3,3> refineFundamentalMatrix(std::vector<Match> & bestRansacMatches){
  FMatrix<float,3,3> F(0.0);

  int maxNbInliers = bestRansacMatches.size();
  
  FMatrix<float,3,3> normalizationMatrix(0.0);
  normalizationMatrix(0,0) = 0.001;
  normalizationMatrix(1,1) = 0.001;
  normalizationMatrix(2,2) = 1;
  
  Matrix<float> rA(maxNbInliers,9);
  Matrix<float> Ura(maxNbInliers,maxNbInliers);
  Vector<float> Sra(maxNbInliers);
  Matrix<float> Vtra(9,9);

  Vector<float> fVector(9);

  FMatrix<float,3,3> reformedF;
  FMatrix<float,3,3> Urf(0.0), Vtrf(0.0);
  FVector<float,3> Srf(0.0);

  
  //Compute A matrix
  for(int i = 0; i < maxNbInliers; i++){
    FVector<float,3> pl(0),pr(0);
    Match m = bestRansacMatches.back();
    bestRansacMatches.pop_back();
    pl[0] = m.x1;
    pl[1] = m.y1;
    pl[2] = 1;
    
    pr[0] = m.x2;
    pr[1] = m.y2;
    pr[2] = 1;
    
    pl = normalizationMatrix*pl;
    pr = normalizationMatrix*pr;
    
    
    rA(i,0) = pl[0]*pr[0] ;
    rA(i,1) = pl[0]*pr[1] ;
    rA(i,2) = pl[0] ;
    rA(i,3) = pl[1]*pr[0] ;
    rA(i,4) = pl[1]*pr[1] ;
    rA(i,5) = pl[1] ;
    rA(i,6) = pr[0] ;
    rA(i,7) = pr[1] ;
    rA(i,8) = 1 ;
  }    
  

  svd(rA,Ura,Sra,Vtra);
  
  fVector = Vtra.getRow(8);
  
  for(int i = 0; i < 9;i++){
    reformedF(i/3,i%3) = fVector[i];
  }
  
  
  svd(reformedF, Urf,Srf,Vtrf);

  //Ensure rank of F matrix
  Srf[2] = 0;
  
  //Form the F matrix back
  F = Urf*(Diagonal(Srf)*Vtrf);
  

  //Unnormalize F matrix
  F = transpose(normalizationMatrix)*F*normalizationMatrix;
  

  //Ensure that the norm = 1 
  F = F/norm(F);

  return F;
}

FMatrix<float,3,3> computeF(const std::vector<Match>& matches) {
    FMatrix<float,3,3> F;
    vector<Match> bestRansacMatches;

    //Compute the inliers of the best fundamental matrix estimate
    bestRansacMatches = computeRansacInliers(matches);
    
    //Refine the fundamental matrix estimate using the inliers
    F = refineFundamentalMatrix(bestRansacMatches);
    
    return F;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const FMatrix<float,3,3>& F) {
    while(true) {
        int x,y;
        if(getMouse(x,y) == 3)
            break;
	if(x<I1.width()){
	  FVector<float,3> pr(0.0), pl(0.0);
	  FMatrix<float,3,3> Ft = transpose(F);
	  pl[0] = x;
	  pl[1] = y;
	  pl[2] = 1; 

	  //compute the right epipolar line equation associated to the point pl
	  pr = Ft*pl;


	  //Compute coordinates for two specific points of the line (x=0;x=Width)
	  FVector<float,3> epl1(0),epl2(0);
	  epl1[0] = 0;
	  epl1[1] = (-pr[2] - pr[0]*epl1[0])/pr[1];
	  
	  epl2[0] = I1.width();
	  epl2[1] = (-pr[2] - pr[0]*epl2[0])/pr[1];

	  epl1[0] = epl1[0] + I1.width();
	  epl2[0] = epl2[0] + I1.width();

	  drawLine(epl1[0],epl1[1],epl2[0],epl2[1],RED,2);
	  
	}
	else{
	  FVector<float,3> pr(0.0), pl(0.0);
	  pr[0] = x-I1.width();
	  pr[1] = y;
	  pr[2] = 1; 

	  pr = pr;
	  

	  //Compute the left epipolar line equation associated the the point pr
	  pl = F*pr;

	  //Compute coordinates for two specific points of the line (x=0;x=Width)
	  FVector<float,3> epl1(0.0), epl2(0.0);

	  epl1[0] = 0;
	  epl1[1] = (-pl[2] - pl[0]*epl1[0])/pl[1];
	  epl2[0] = I1.width();
	  epl2[1] = (-pl[2] - pl[0]*epl2[0])/pl[1];

	  drawLine(epl1[0],epl1[1],epl2[0],epl2[1],RED,2);
	}
    }
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));

    const char* s1 = argc>1? argv[1]: srcPath("im1.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("im2.jpg");

    // Load and display images
    Image<Color,2> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load images" << endl;
        return 1;
    }
    int w = I1.width();
    openWindow(2*w, I1.height());
    display(I1,0,0);
    display(I2,w,0);

    std::vector<Match> matches;
    algoSIFT(I1, I2, matches);
    cout << " matches: " << matches.size() << endl;
    click();
    
    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);

    FMatrix<float,3,3> F = computeF(matches);

    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}
