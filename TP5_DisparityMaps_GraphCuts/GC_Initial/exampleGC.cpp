#include <Imagine/Images.h>
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <Imagine/Common.h>
#include <iostream>
#include <fstream>

#include "maxflow/graph.h"

using namespace std;
using namespace Imagine;


void binRestoration(){
  Image<byte> I;
  load(I,srcPath("lena.jpg"));
  openWindow(I.width(), I.height());
  for(int i = 0; i < 300;i++){
    for(int j = 0; j < 300;j ++){
      if(I(i,j) > 128){
	I(i,j) = 255;
      }
      else{
	I(i,j) = 0;
      }
    }
  }
  display(I);
  
  //Graph allocation
  Graph<int,int,int> g(900000,2*900000 /*sing + source*/ + 300*300*4 /*interconnections*/ - (300+300+298+298)) /*borders*/;
  g.add_node(900000);
  
  int k = 10;
  
  //Graph penalty weights
  for(int i = 0; i < 300;i++){
    for(int j = 0; j < 300; j++){
      if(I(i,j)==255){
	g.add_tweights(i*300+j,k,0);
      }
      else{
	g.add_tweights(i*300+j,0,k);
      }
    }
  }
  
  int lambda = 10;
  //Graph regularisation weigths
  for(int i = 0; i < 300; i++){
    for(int j = 0; j < 300; j++){
      if(j!=299)
	g.add_edge(i*300+j, i*300+j+1,lambda,lambda); // right;
      /*if(j!=0)
	g.add_edge(i*300+j, i*300+j-1,lambda,lambda); // left;*/
      if(i!=299)
	g.add_edge(i*300+j, i*300+j+300,lambda,lambda); // bottom;
      /*if(i!=0)
	g.add_edge(i*300+j, i*300+j-300,lambda,lambda); // top;*/
      
    }
  }
  
  int flow = g.maxflow();
  
  Image<byte> I2(300,300);
  
  for(int i = 0; i < 300; i++){
    for(int j = 0; j < 300; j++){
      if(g.what_segment(i*300+j)==Graph<int,int,int>::SOURCE){
	I2(i,j) = 255;
      }
      else{
	I2(i,j) = 0;
      }
    }
  }
  cout << "Flow : " << flow << endl;
  
  Window w2 = openWindow(I2.width(),I2.height());
  setActiveWindow(w2);
  display(I2);
  
}


void segmentImages()
{
  Image<Color> I1;
  load(I1,srcPath("black_horse.jpg"));

  vector<IntPoint2> foregroundPoints;
  vector<IntPoint2> backgroundPoints;
  bool finishedClickingForeground = false;
  bool finishedClickingBackground = false;
  
  Window w = openWindow(I1.width(), I1.height());
  setActiveWindow(w);
  display(I1);
  
  //Wait for clicks for foreground images
  cout << "Waiting for clicks in the foreground" << endl;
  IntPoint2 p;
  while(!finishedClickingForeground){
    Window win;
    int subwin;
    int button = anyGetMouse(p,win,subwin);
    if(button==3){
      finishedClickingForeground = true;
      continue;
    }
    foregroundPoints.push_back(p);
    fillCircle(p,3,RED);
    
  }

  //Wait for clicks for background images  
  cout << "Waiting for clicks in the background" << endl;
  while(!finishedClickingBackground){
    Window win;
    int subwin;
    int button = anyGetMouse(p,win,subwin);
    if(button==3){
      finishedClickingBackground = true;
      continue;
    }
    backgroundPoints.push_back(p);
    fillCircle(p,3,BLUE);
  }
  
  //Compute color of object
  float fgRed = 0,fgGreen = 0,fgBlue = 0;
  float sigmaFgRed = 0, sigmaFgGreen = 0, sigmaFgBlue = 0;
  
  float bgRed = 0,bgGreen = 0,bgBlue = 0;
  float sigmaBgRed = 0, sigmaBgGreen = 0, sigmaBgBlue = 0;
  
  //mean = 1/N sum(xi)
  //std = sqrt(1/N*sum((x-mx)2))

  float foregroundSize = foregroundPoints.size();
  float backgroundSize = backgroundPoints.size();
  
  int count = 0;
  for(int i = 0; i < foregroundSize; i++){
    cout << foregroundPoints[i] << endl;
    for(int k = -1; k <=1; k++){
      for(int l = -1; l <= 1; l++){
	if(foregroundPoints[i].x()+k >= 0 && foregroundPoints[i].x()+k<I1.width() &&
	   foregroundPoints[i].y()+l >= 0 && foregroundPoints[i].y()+l<I1.height()){
	  Color c = I1(foregroundPoints[i].x()+k, foregroundPoints[i].y()+l);
	  cout << c << endl;
	  fgRed += c[0];
	  fgGreen += c[1];
	  fgBlue += c[2];
	  cout << fgBlue << endl;
	  count++;
	}
      }
    }
  }

  fgRed /= count;
  fgGreen /= count;
  fgBlue /= count;

  count = 0;
  for(int i = 0; i < foregroundSize; i++){
    for(int k = -1; k <= 1; k++){
      for(int l = -1; l <= 1; l++){
	if(foregroundPoints[i].x()+k >= 0 && foregroundPoints[i].x()+k<I1.width() &&
	   foregroundPoints[i].y()+l >= 0 && foregroundPoints[i].y()+l<I1.height()){
	  Color c = I1(foregroundPoints[i].x()+k, foregroundPoints[i].y()+l);
	  sigmaFgRed += (c[0] - fgRed)*(c[0] - fgRed);
	  sigmaFgGreen += (c[1] - fgGreen)*(c[1] - fgGreen);
	  sigmaFgBlue += (c[2] - fgBlue)*(c[2] - fgBlue);
	  count++;
	}
      }
    } 
  }

  sigmaFgRed = sqrt(sigmaFgRed/count);
  sigmaFgGreen = sqrt(sigmaFgGreen/count);
  sigmaFgBlue = sqrt(sigmaFgBlue/count);

  //
  cout << endl;
  count = 0;
  for(int i = 0; i < backgroundSize; i++){
    cout << backgroundPoints[i] << endl;
    for(int k = -1; k <=1; k++){
      for(int l = -1; l <= 1; l++){
	if(backgroundPoints[i].x()+k >= 0 && backgroundPoints[i].x()+k<I1.width() &&
	   backgroundPoints[i].y()+l >= 0 && backgroundPoints[i].y()+l<I1.height()){
	  Color c = I1(backgroundPoints[i].x()+k, backgroundPoints[i].y()+l);
	  cout << c << endl;
	  bgRed += c[0];
	  bgGreen += c[1];
	  bgBlue += c[2];
	  count++;
	}
      }
    }
  }

  bgRed /= count;
  bgGreen /= count;
  bgBlue /= count;
  

  count = 0;
  for(int i = 0; i < backgroundSize; i++){
    for(int k = -1; k <= 1; k++){
      for(int l = -1; l <= 1; l++){
	if(backgroundPoints[i].x()+k >= 0 && backgroundPoints[i].x()+k<I1.width() &&
	   backgroundPoints[i].y()+l >= 0 && backgroundPoints[i].y()+l<I1.height()){
	  Color c = I1(backgroundPoints[i].x()+k, backgroundPoints[i].y()+l);
	  sigmaBgRed += (c[0] - bgRed)*(c[0] - bgRed);
	  sigmaBgGreen += (c[1] - bgGreen)*(c[1] - bgGreen);
	  sigmaBgBlue += (c[2] - bgBlue)*(c[2] - bgBlue);
	  count++;
	}
      }
    } 
  }

  sigmaBgRed = sqrt(sigmaBgRed/count);
  sigmaBgGreen = sqrt(sigmaBgGreen/count);
  sigmaBgBlue = sqrt(sigmaBgBlue/count);

  

  cout << "Red : (" <<fgRed << "," << sigmaFgRed<<")" << endl;
  cout << "Green : (" << fgGreen << "," << sigmaFgGreen << ")" << endl;
  cout << "Blue : (" << fgBlue << "," << sigmaFgBlue << ")" << endl;

  cout << "Red : (" <<bgRed << "," << sigmaBgRed<<")" << endl;
  cout << "Green : (" << bgGreen << "," << sigmaBgGreen << ")" << endl;
  cout << "Blue : (" << bgBlue << "," << sigmaBgBlue << ")" << endl;
  
  Matrix<float> muFG(3,1);
  muFG(0,0) = fgRed;
  muFG(1,0) = fgGreen;
  muFG(2,0) = fgBlue;

  Matrix<float> muBG(3,1);
  muBG(0,0) = bgRed;
  muBG(1,0) = bgGreen;
  muBG(2,0) = bgBlue;

  Matrix<float> sFG(3,3);
  sFG(0,0) = sigmaFgRed;
  sFG(0,1) = 0;
  sFG(0,2) = 0;
  sFG(1,0) = 0;
  sFG(1,1) = sigmaFgGreen;
  sFG(1,2) = 0;
  sFG(2,0) = 0;
  sFG(2,1) = 0;
  sFG(2,2) = sigmaFgBlue;


  
  Matrix<float> sBG(3,3);
  sBG(0,0) = sigmaBgRed;
  sBG(0,1) = 0;
  sBG(0,2) = 0;
  sBG(1,0) = 0;
  sBG(1,1) = sigmaBgGreen;
  sBG(1,2) = 0;
  sBG(2,0) = 0;
  sBG(2,1) = 0;
  sBG(2,2) = sigmaBgBlue;
  
  
  
  //Create graph
  //Source : foreground
  //Sink : background
  //Size : nb pixels dans l'image
  //terme d'attache aux données : distance à la gaussienne
  //terme de régularisation : lambda
  
  int nbNodes = I1.width()*I1.height();
  int nbEdges = nbNodes*6;
  float lambda = 2.;
  float sigma = 5.;

  Matrix<float> sigmas(3,3);
  
  sigmas(0,0) = sigma;
  sigmas(0,1) = 0;
  sigmas(0,2) = 0;
  sigmas(1,0) = 0;
  sigmas(1,1) = sigma;
  sigmas(1,2) = 0;
  sigmas(2,0) = 0;
  sigmas(2,1) = 0;
  sigmas(2,2) = sigma;
  

  cout << "Creating graph" << endl;
  Graph<float,float,float> g(nbNodes,nbEdges);
  g.add_node(nbNodes);

  cout << "Asigning weights" << endl;
  for(int j = 0; j < I1.height(); j++){
    for(int i = 0; i < I1.width(); i++){
      Color c = I1(i,j);
      
      float probaForeground = 1/(sqrt(pow(2*3.14159256,3))*sqrt(norm(sFG)));
      float probaBackground = 1/(sqrt(pow(2*3.14159256,3))*sqrt(norm(sBG)));
      
      Matrix<float> x(3,1);
      x(0,0) = c[0];
      x(1,0) = c[1];
      x(2,0) = c[2];
        
      probaForeground *= exp(-0.5*(transpose(x-muFG)*inverse(sFG)*(x-muFG))[0]);
      probaBackground *= exp(-0.5*(transpose(x-muBG)*inverse(sBG)*(x-muBG))[0]);
      

      g.add_tweights(j*I1.width()+i,-log(probaForeground),-log(probaBackground));

      if(j<I1.height()-1){//On peut ajouter une edge vers le bas
       
	float B = 0;
	Color c1 = I1(i,j);
	Color c2 = I1(i,j+1);
	
	Matrix<float> col1(3,1),col2(3,1);
	
	col1(0,0) = c1[0];
	col1(1,0) = c1[1];
	col1(2,0) = c1[2];
	
	col2(0,0) = c2[0];
	col2(1,0) = c2[1];
	col2(2,0) = c2[2];
	


	B = exp(-0.5*(transpose(col1-col2)*inverse(sigmas)*(col1-col2))[0]);
	
	
	g.add_edge(j*I1.width()+i ,(j+1)*I1.width()+i,lambda,lambda);
      }
      
      if(i<I1.width()-1){//On peut ajouter une edge vers la droite
	g.add_edge(j*I1.width()+i, j*I1.width()+i+1,lambda,lambda);
      }
    }
  }
  
  
  cout << "Computing max flow" << endl;
  int maxFlow = g.maxflow();
  cout << "Max flow : " << maxFlow << endl;
  //maxflow graph
  
  Image<byte> output(I1.width(), I1.height());
  for(int j = 0; j < I1.height(); j++){
    for(int i = 0; i < I1.width(); i++){
      if(g.what_segment(j*I1.width()+i)==Graph<float,float,float>::SOURCE){
	output(i,j) = 255;
      }
      else{
	output(i,j) = 0;
      }
    }
  }

  //output = blur(output,3.);

  display(output);

  
  //create result image
}


int main() {
	binRestoration();
	segmentImages();
	endGraphics();
	return 0;
}
