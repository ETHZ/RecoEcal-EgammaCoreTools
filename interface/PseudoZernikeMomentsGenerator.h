#ifndef RecoEcal_EgammaCoreTools_PseudoZernikeMomentsGenerator_h
#define RecoEcal_EgammaCoreTools_PseudoZernikeMomentsGenerator_h

#define PI 3.14159265

#include <map>

class PseudoZernikeMomentsGenerator
{
  public:
    enum energyScalingTypes { NONE, ENERGYSUM, SEEDCRYSTAL }; 
    
    PseudoZernikeMomentsGenerator(int radius, std::map<std::pair<int,int>,double> digitizedEnergyMap, bool useLogWeightedMoments, double param_W0, bool useRotationInvariantMap, energyScalingTypes energyScalingType);
    ~PseudoZernikeMomentsGenerator();
    
    double calculateMoment(int whichMoment); 
    std::map<std::pair<int,int>,double> getDigitizedEnergyMap() { return digitizedEnergyMap_; }
   
  private:

    void findEnergySumForScaleInvariance(energyScalingTypes energyScalingType);  
    void makeEnergyMapTransformationInvariant();
    void rotateEnergyMap(bool rotateCounterClockwise);
    void mirrorEnergyMap(bool mirrorVertically);        
    
    double calculateRDependence(int n, int m, double normalizedRadius); 
    double calculateThetaDependence(int m, int p, double theta);
    
    inline double findNormalizedRadius(int x, int y); 
    inline double findTheta(int x, int y);
           double findFactorial(int number); 
    
    int radius_; 
    double energySum_;
    std::map<std::pair<int,int>,double> digitizedEnergyMap_;
    bool useLogWeightedMoments_;
    double param_W0_;
};

#endif
