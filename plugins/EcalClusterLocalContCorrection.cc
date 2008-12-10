#include "RecoEcal/EgammaCoreTools/plugins/EcalClusterLocalContCorrection.h"


float EcalClusterLocalContCorrection::getValue( const reco::BasicCluster & basicCluster, const EcalRecHitCollection & recHit) const
{
        checkInit();
        // private member params_ = EcalClusterLocalContCorrectionParameters
        // (see in CondFormats/EcalObjects/interface)
        EcalClusterLocalContCorrParameters::const_iterator it;
        std::cout << "[[EcalClusterLocalContCorrection::getValue]] " 
                << params_->size() << " parameters:";
        for ( it = params_->begin(); it != params_->end(); ++it ) {
                std::cout << " " << *it;
        }
        std::cout << "\n";
        return 1;
}

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionFactory.h"
DEFINE_EDM_PLUGIN( EcalClusterFunctionFactory, EcalClusterLocalContCorrection, "EcalClusterLocalContCorrection");
