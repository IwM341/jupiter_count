#include "functions.hpp"
#include <cmath>
#include <fstream>

int main(void){
	double mk = 0.7;
	double mp = 0.938;
	
	auto PhiFactor = [mp,mk](double dv2)->double{
			return dv2*dv2*(mk*mk/mp*mp);
		};
	
	const auto& BM = BodyModel(2.03e-4,"../DMFramework/celstial_models/jupiter_model.dat");   
	
	SupressFactor(mk,ELASTIC,PhiFactor,BM,"H",100000,U0*1.1,U0);
	SupressFactor(mk,IONIZATION,PhiFactor,BM,"H",100000,U0*1.1,U0);
	SupressFactor(mk,MIGDAL,PhiFactor,BM,"H",100000,U0*1.1,U0);
	
	std::vector<double> Vmk = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8};
	size_t N = Vmk.size();
	
	std::vector<double> Fel(N);
	std::vector<double> Fmgd(N);
	std::vector<double> Fion(N);
	
	#pragma omp parallel for
	for(size_t i=0;i<N;++i){
		Fel[i] = SupressFactor(Vmk[i],ELASTIC,PhiFactor,BM,"H",1000000,U0*1.1,U0);
		Fion[i] = SupressFactor(Vmk[i],IONIZATION,PhiFactor,BM,"H",1000000,U0*1.1,U0);
		Fmgd[i] = SupressFactor(Vmk[i],MIGDAL,PhiFactor,BM,"H",1000000,U0*1.1,U0);
	}
	
	std::ofstream ofsEl("elastic.dat");
	std::ofstream ofsIon("ion.dat");
	std::ofstream ofsMgd("migdal.dat");
	
	
	ofsEl << "mk\tF\n"<< Function::GridFunction1(Vmk,Fel) << std::endl;
	ofsIon << "mk\tF\n"<< Function::GridFunction1(Vmk,Fion) << std::endl;
	ofsMgd << "mk\tF\n"<< Function::GridFunction1(Vmk,Fmgd) << std::endl;
	
	return 0;
	
}
