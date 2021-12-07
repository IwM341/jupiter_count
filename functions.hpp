#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <utils>
#include <random>
#include <cmath>

#define U0 0.7667e-3
#define JVESC 2.03e-4
#define alpha 0.0073

template <class Generator> 
inline double RandomCos(Generator G){
	return G()*2.0-1;
}
template <class Generator> 
inline double RandomPhi(Generator G){
	return G()*2*M_PI;
}

template <class Generator> 
inline MC::MCResult<vec3> Boltsman(Generator G, double Vdisp,double Vmin = 0){
	double V = sqrt(Vmin*Vmin-2*Vdisp*Vdisp*log(G()));
	return MC::MCResult<vec3>(vec3::PolarCos(V,RandomCos(G),0),
				exp(-Vmin*Vmin/(2*Vdisp*Vdisp))*sqrt(2/M_PI)*V/Vdisp);
}

template <class Generator> 
inline MC::MCResult<vec3> Velocity(Generator G,double VescTmp,
						double Vdisp = 1.1*U0,double mU0 = U0){
	double u_nd = sqrt(-2*log(G()) );
	double U = u_nd*Vdisp;
	double u0_nd = mU0/Vdisp;
	
	double V = sqrt(U*U+VescTmp*VescTmp);
	
	double rd = V/Vdisp*exp(0.5*u0_nd*(2*u_nd-u0_nd))*
				(1-exp(-2*u0_nd*u_nd))/(u_nd*sqrt(2*M_PI));
	
	return MC::MCResult<vec3>(vec3(V,0,0),rd);
	
}

template <class Generator> 
inline MC::MCResult<vec3> NuOut(Generator G,const vec3& Vcm,const vec3&Nu,
						double Vesc,double mp,double mk,double deltaE = 0){
	double VcmN = Vcm.norm();
	
	vec3 n_v = Vcm/VcmN;
	vec3 n_1(-n_v.z,0,n_v.x);
	vec3 n_y(0,1,0);
	
	double Nu1_squared = 	Nu.quad()-deltaE*2*mp/(mk*(mp+mk));
	if(Nu1_squared<=0.0)
		return MC::MCResult<vec3>(vec3(0,0,0),0);
	
	double Nu1 = sqrt(Nu1_squared);
	
	double cosTh1max = (Vesc*Vesc-Nu1_squared-VcmN*VcmN)/(2*VcmN*Nu1);
	
	if(cosTh1max <= -1)
		return MC::MCResult<vec3>(vec3(0,0,0),0);
	else if(cosTh1max >= 1){
		cosTh1max = 1;
	}
	
	double cosTh1 = (1+cosTh1max)*G()-1;
	
	return MC::MCResult<vec3>(vec3::PolarCos(Nu1,cosTh1,RandomPhi(G)),
								0.5*(1.0+cosTh1max)*Nu1);
	
}

enum ScatteringType{
	ELASTIC,IONIZATION,MIGDAL
};
template <typename FuncType>
double SupressFactor(double mk,ScatteringType ST,
					FuncType PhiFactor, 
					const BodyModel& BM,const std::string &element,
					size_t Nmk,
					double Vdisp,double mU0 = U0){
	
	std::default_random_engine stdGen;
	std::uniform_real_distribution<double> stdUniform(1.0,0.0);
	auto G = [&stdGen,&stdUniform]()->double{return stdUniform(stdGen);};
	
	const double mp = 0.938;
	
	double sum = 0;
	double sum2 = 0;
	
	auto VescR = BM.UniformRadFunc("Vesc");
	
	for(size_t i=0;i<Nmk;++i){
		
		double factor = 1.0;
		
		double r_nd = pow(G(),1.0/3.0);
		
		double Vesc = VescR(r_nd);
		auto VelocityMk = Velocity(G,Vesc,Vdisp,mU0);
		
		vec3 V_wimp = VelocityMk.Result;
		factor *= VelocityMk.RemainDensity;
		
		double n_nd = 1.0;//TODO n_nd as a function of radius
		
		vec3 V1(0,0,0);//TODO: add thermal distribution of nuclei velocity
		
		
		double deltaE = 0;
		double Inelastic_rd = 1.0;
		if(ST == IONIZATION){
			auto Rs = PhiGenerate(G);
			deltaE = deltaEIon(Rs.Result);
			Inelastic_rd = Rs.RemainDensity;
		}
		else if(ST == MIGDAL){
			auto Rs = NGenerate(G);
			deltaE = deltaEMgd(Rs.Result);
			Inelastic_rd = Rs.RemainDensity;
		}
		
		factor*= Inelastic_rd;
		
		vec3 Vcm = (V_wimp*mk + V1*mp)/(mp+mk);
		vec3 Nu = mp/(mp+mk)*(V_wimp-V1);
		
		auto Numk = NuOut(G,Vcm,Nu,Vesc,mp,mk,deltaE);
		
		vec3 Nu1 = Numk.Result;
		factor*=Numk.RemainDensity;
		
		double dv2 = (Nu-Nu1).quad();
		
		if(ST != ELASTIC){
			factor *= dv2*pow(mk/(alpha*mp),2);
		}
		
		factor *= PhiFactor(dv2);
		
		sum += factor/Nmk;
		sum2 += factor*factor/Nmk;
		
		
	}
	/*
	std::cout << "integration competed, sum = " << sum << ", disp = " << 
					sqrt(sum2-sum*sum) << ", precision = " << 
					sqrt( (sum2-sum*sum)/(Nmk-1)) << std::endl;
	*/
	return sum;
	
}





#endif
