#include "RecoEcal/EgammaCoreTools/interface/PseudoZernikeMomentsGenerator.h"
#include <cmath>

PseudoZernikeMomentsGenerator::PseudoZernikeMomentsGenerator(int radius, std::map<std::pair<int,int>,double> digitizedEnergyMap, bool useLogWeightedMoments, double param_W0, bool useRotationInvariantMap, energyScalingTypes energyScalingType) : 
  radius_(radius), digitizedEnergyMap_(digitizedEnergyMap), useLogWeightedMoments_(useLogWeightedMoments), param_W0_(param_W0)
{
  findEnergySumForScaleInvariance(energyScalingType);
  if(useRotationInvariantMap) makeEnergyMapTransformationInvariant();
}

PseudoZernikeMomentsGenerator::~PseudoZernikeMomentsGenerator() {}

void PseudoZernikeMomentsGenerator::findEnergySumForScaleInvariance(energyScalingTypes energyScalingType)
{
  if(energyScalingType == NONE) energySum_ = 1; 
  else
  {
    energySum_ = 0; 

    if(energyScalingType == ENERGYSUM)
    {
      for(int x = -1*radius_; x <= radius_; x++)  
        for(int y = -1*radius_; y <= radius_; y++)
          energySum_ += (useLogWeightedMoments_) ? std::max(0., param_W0_ + log(digitizedEnergyMap_.find(std::pair<int,int>(x,y))->second)) : digitizedEnergyMap_.find(std::pair<int,int>(x,y))->second;
    }
    else if(energyScalingType == SEEDCRYSTAL)
      energySum_ = (useLogWeightedMoments_) ? std::max(0., param_W0_ + log(digitizedEnergyMap_.find(std::pair<int,int>(0,0))->second)) : digitizedEnergyMap_.find(std::pair<int,int>(0,0))->second;
  }
}

void PseudoZernikeMomentsGenerator::makeEnergyMapTransformationInvariant()
{
  double northEnegry = 0;
  double eastEnegry = 0;
  double southEnegry = 0;
  double westEnegry = 0;

  for(int i = 1; i < radius_; i++)
  {
    northEnegry += digitizedEnergyMap_.find(std::pair<int,int>(0,-1*i))->second;
    eastEnegry  += digitizedEnergyMap_.find(std::pair<int,int>(1*i,0))->second;
    southEnegry += digitizedEnergyMap_.find(std::pair<int,int>(0,1*i))->second;
    westEnegry  += digitizedEnergyMap_.find(std::pair<int,int>(-1*i,0))->second;            
  }

  enum {North, East, South, West};
  std::map<double,int> cardinalEnergyMap; 
  cardinalEnergyMap.insert(std::make_pair(northEnegry,(int)North));
  cardinalEnergyMap.insert(std::make_pair(eastEnegry,(int)East));
  cardinalEnergyMap.insert(std::make_pair(southEnegry,(int)South));
  cardinalEnergyMap.insert(std::make_pair(westEnegry,(int)West));    
  
  std::map<double,int>::reverse_iterator cardinalEnergyMapIterator = cardinalEnergyMap.rbegin();
  int maxCardinalEnergyDirection = cardinalEnergyMapIterator->second;
  cardinalEnergyMapIterator++;
  int secondCardinalEnergyDirection = cardinalEnergyMapIterator->second;
  
  if(maxCardinalEnergyDirection == East && secondCardinalEnergyDirection == North)       {rotateEnergyMap(true);}
  else if(maxCardinalEnergyDirection == South && secondCardinalEnergyDirection == East)  {rotateEnergyMap(true); rotateEnergyMap(true);}  
  else if(maxCardinalEnergyDirection == South && secondCardinalEnergyDirection == West)  {mirrorEnergyMap(true);}    
  else if(maxCardinalEnergyDirection == West && secondCardinalEnergyDirection == South)  {rotateEnergyMap(false);}    
  else if(maxCardinalEnergyDirection == West && secondCardinalEnergyDirection == North)  {rotateEnergyMap(false);}    
  else if(maxCardinalEnergyDirection == North && secondCardinalEnergyDirection == East)  {mirrorEnergyMap(false);}    
  else if(maxCardinalEnergyDirection == East && secondCardinalEnergyDirection == South)  {rotateEnergyMap(true);}      
  else if(maxCardinalEnergyDirection == South && secondCardinalEnergyDirection == North) {mirrorEnergyMap(true);}    
  else if(maxCardinalEnergyDirection == East && secondCardinalEnergyDirection == West)   {rotateEnergyMap(true);}    
  else if(maxCardinalEnergyDirection == West && secondCardinalEnergyDirection == East)   {rotateEnergyMap(false);}
  
  if(digitizedEnergyMap_.find(std::pair<int,int>(1,0))->second > digitizedEnergyMap_.find(std::pair<int,int>(-1,0))->second) mirrorEnergyMap(false); 

  cardinalEnergyMap.clear(); 
  cardinalEnergyMap.insert(std::make_pair(digitizedEnergyMap_.find(std::pair<int,int>(0,-1))->second,(int)North));
  cardinalEnergyMap.insert(std::make_pair(digitizedEnergyMap_.find(std::pair<int,int>(1,0))->second,(int)East));
  cardinalEnergyMap.insert(std::make_pair(digitizedEnergyMap_.find(std::pair<int,int>(0,1))->second,(int)South));
  cardinalEnergyMap.insert(std::make_pair(digitizedEnergyMap_.find(std::pair<int,int>(-1,0))->second,(int)West));    
}

void PseudoZernikeMomentsGenerator::rotateEnergyMap(bool rotateCounterClockwise)
{
  std::map<std::pair<int,int>,double> copyOfDigitizedEnergyMap(digitizedEnergyMap_);
  digitizedEnergyMap_.clear();
  
  for(int x = -1*radius_; x <= radius_; x++) 
    for(int y = -1*radius_; y <= radius_; y++)
      if(rotateCounterClockwise) digitizedEnergyMap_.insert(std::make_pair(std::pair<int,int>(x,y),copyOfDigitizedEnergyMap.find(std::pair<int,int>(-1*y,x))->second));
      else digitizedEnergyMap_.insert(std::make_pair(std::pair<int,int>(x,y),copyOfDigitizedEnergyMap.find(std::pair<int,int>(y,-1*x))->second));        
}

void PseudoZernikeMomentsGenerator::mirrorEnergyMap(bool mirrorVertically)
{
  std::map<std::pair<int,int>,double> copyOfDigitizedEnergyMap(digitizedEnergyMap_);
  digitizedEnergyMap_.clear();
  
  for(int x = -1*radius_; x <= radius_; x++) 
    for(int y = -1*radius_; y <= radius_; y++)
      if(mirrorVertically) digitizedEnergyMap_.insert(std::make_pair(std::pair<int,int>(x,y),copyOfDigitizedEnergyMap.find(std::pair<int,int>(x,-1*y))->second));
      else  digitizedEnergyMap_.insert(std::make_pair(std::pair<int,int>(x,y),copyOfDigitizedEnergyMap.find(std::pair<int,int>(-1*x,y))->second));
}
  
double PseudoZernikeMomentsGenerator::calculateMoment(int whichMoment)
{
  int n = (int)std::floor(std::sqrt(whichMoment));
  int m = (int)(n - std::floor((whichMoment-std::pow((double)n,2))/2.));
  int p = (m!=0) ? (int)(whichMoment-std::pow((double)n,2))%2 : -1;
  
  double pseudoZernikeMoment = 0.0;

  for(int x = -1*radius_; x <= radius_; x++)
    for(int y = -1*radius_; y <= radius_; y++)
    {
      double theta = findTheta(x,y);
      double normalizedRadius = findNormalizedRadius(x,y);
      
      pseudoZernikeMoment += (useLogWeightedMoments_) ? 
        ((double)(m+1.0)/PI) * std::max(0., param_W0_ + log(digitizedEnergyMap_.find(std::pair<int,int>(x,y))->second)) * calculateRDependence(n,m,normalizedRadius) * calculateThetaDependence(n,p,theta) :
        ((double)(m+1.0)/PI) * digitizedEnergyMap_.find(std::pair<int,int>(x,y))->second * calculateRDependence(n,m,normalizedRadius) * calculateThetaDependence(n,p,theta);           
    }

  return pseudoZernikeMoment/energySum_;
}

double PseudoZernikeMomentsGenerator::findNormalizedRadius(int x, int y)
{
  return std::sqrt(std::pow((double)x,2)+std::pow((double)y,2))/(radius_);
}

double PseudoZernikeMomentsGenerator::findTheta(int x, int y)
{ 
  return std::atan2((double)y,(double)x); 
}

double PseudoZernikeMomentsGenerator::findFactorial(int number) 
{
	if(number <= 1) return 1;
	return number * findFactorial(number - 1);
}

double PseudoZernikeMomentsGenerator::calculateRDependence(int n, int m, double normalizedRadius)
{
  double calculatedRDependence = 0.0;   

  for(int s = 0; s <= n-m; s++)
    calculatedRDependence +=  std::pow(-1.,s)*findFactorial(2*n-m-s)
                              /(findFactorial(s)*findFactorial(n-s)*findFactorial(n-m-s))
                              *std::pow(normalizedRadius,2*(n-s)-m);

  return calculatedRDependence;
}

double PseudoZernikeMomentsGenerator::calculateThetaDependence(int m, int p, double theta)
{
  switch(p)
  {
    case 0: return std::cos(m*theta); break;
    case 1: return std::sin(m*theta); break;
    default: return 1;
  }
}

 
