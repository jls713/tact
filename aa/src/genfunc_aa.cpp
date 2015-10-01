// ============================================================================
// Action calculation using generating function from orbit integration
// ============================================================================

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "GSLInterface/GSLInterface.h"
#include "utils.h"
#include "potential.h"
#include "analytic_aa.h"
#include "stackel_aa.h"
#include "genfunc_aa.h"
#include "debug.h"
#include "orbit.h"

static int sign(double a){
	if(a>=0.) return 1;
	else return -1;
}

static void orient_orbit(VecDoub& x, const std::vector<int> &loop){
    // given which sign of the angular momentum is conserved, orients the orbit
    // appropriately
	if(loop[0]){ // L_x sign preserved
		double tmp = x[0];x[0]=-x[2];x[2]=tmp;
        tmp = x[3];x[3]=-x[5];x[5]=tmp;
	}
	if(loop[1]){ // L_y sign preserved
		double tmp = x[1];x[1]=-x[2];x[2]=tmp;
        tmp = x[4];x[4]=-x[5];x[5]=tmp;
	}
}

std::vector<int> Actions_Genfunc::angular_momentum(const VecDoub &x){
	VecDoub xx = {x[0],x[1],x[2]}, vv = {x[3],x[4],x[5]};
	VecDoub ll = cross_product<double>(xx,vv);
	return {sign(ll[0]),sign(ll[1]),sign(ll[2])};
}

std::vector<int> Actions_Genfunc::loop(const std::vector<VecDoub> &orbit_samples){
    // from a series of phase-space points assesses whether signs of angular
    // momentum are conserved -- returns \pm 1 depending on sign of conserved
    // ang mom or zero otherwise
	std::vector<int> l(3,true);
	std::vector<int> l0=angular_momentum(orbit_samples[0]);
    std::vector<int> result = l0;
	for(auto it = std::begin(orbit_samples)+1;
	    	 it!=std::end(orbit_samples);
	    	 ++it)
	{
		l = angular_momentum(*it);
		for(int i=0;i<3;i++) if(l[i]!=l0[i]) result[i] = 0;
		if(result[0]==0 and result[1]==0 and result[2]==0) break;
	}
	return result;
}

// Interface with LAPACK
extern "C" {
// extern void dsptrf_(char *,int *,double *,int *,int*);
// extern void dsptrs_(char *,int *,int *,double *,int*,double*,int*,int*);
extern void dspsv_(char *,int *,int *,double *,int*,double*,int*,int*);
// extern void dsysv_(char *,int *,int *,double *,int*,double*,int*,double*,int*,int*);
};

struct isoparams_st{
    double r1, r2, phi1, phi2;
    isoparams_st(double r1, double r2, double phi1, double phi2)
        :r1(r1),r2(r2),phi1(phi1),phi2(phi2){}
};


double iso_params(double b, void *params){
    /* for finding roots of p_tau^2*2.0*(tau+Alpha)  */
    isoparams_st *RS = (isoparams_st *) params;
    double r2b = sqrt(RS->r2*RS->r2+b*b), r1b = sqrt(RS->r1*RS->r1+b*b);
    return (b+r2b)*(b+r2b)*r2b/((b+r1b)*(b+r1b)*r1b)-RS->phi1*RS->r2/RS->phi2/RS->r1;
}

struct orb{
    const std::vector<VecDoub> &r;
    orb(const std::vector<VecDoub> &r)
        :r(r){}
};


static double fitbox(const gsl_vector *v, void *params){
    orb *O = (orb *) params;
    VecDoub Om(3,0);
    for(unsigned i =0;i<3;++i)
        Om[i] = gsl_vector_get(v,i);
    HarmonicOscillator H0(Om);
    double H = 0., Hm=0.,Hs=0.;
    for(auto i:O->r){
        H=H0.H(i);
        Hm+=H;Hs+=H*H;
    }
    H = (double)(O->r.size());
    Hm/=H;
    return sqrt(Hs/H-Hm*Hm);
}

static double fitiso(const gsl_vector *v, void *params){
    orb *O = (orb *) params;
    VecDoub Om(2,0);
    for(unsigned i =0;i<2;++i)
        Om[i] = gsl_vector_get(v,i);
    Isochrone H0(sqrt(Om[0]),sqrt(Om[1]));
    double H = 0., Hm=0.,Hs=0.;
    for(auto i:O->r){
        H=H0.H(i);
        Hm+=H;Hs+=H*H;
    }
    H = (double)(O->r.size());
    Hm/=H;
    return sqrt(Hs/H-Hm*Hm);
}

VecDoub Actions_Genfunc::find_box_params_minvar(const std::vector<VecDoub> &orbit_samples){
    orb O(orbit_samples);
    minimiser min(&fitbox,{1.,1.,1.},{0.3,0.3,0.3},1e-3,&O);
    VecDoub r(3,0.);
    min.minimise(&r[0],100,false);
    return r;
}

VecDoub Actions_Genfunc::find_isochrone_params_minvar(const std::vector<VecDoub> &orbit_samples){
    orb O(orbit_samples);
    minimiser min(&fitiso,{1.,1.},{0.3,0.3},1e-3,&O);
    VecDoub r(2,0.);
    min.minimise(&r[0],100,false);
    r[0]=sqrt(r[0]);r[1]=sqrt(r[1]);
    return r;
}

VecDoub Actions_Genfunc::find_box_params(const std::vector<VecDoub> &orbit_samples){
    VecDoub xmax(3,0);
    for(auto i:orbit_samples){
        for(int j=0;j<3;j++) if(i[j]>xmax[j]) xmax[j]=i[j];
    }
    VecDoub x(3,0); double F;
    for(int i=0;i<3;i++){
        x[i]=xmax[i];F = -TargetPot->Forces(x)[i];
        xmax[i] = sqrt(F/xmax[i]);
        x[i]=0.;
    }
    return xmax;
}

VecDoub Actions_Genfunc::find_isochrone_params(const std::vector<VecDoub> &orbit_samples){
    // Match inner and outer potential gradients to isochrone
    double r, rmin=1e8, rmax=0.;
    VecDoub in(3,0), out(3,0);
    for(auto i:orbit_samples){
        r=i[0]*i[0]+i[1]*i[1]+i[2]*i[2];
        if(r<rmin){rmin=r;for(int j=0;j<3;j++) in[j]=i[j];}
        if(r>rmax){rmax=r;for(int j=0;j<3;j++) out[j]=i[j];}
    }
    // Fudge factors here an attempt to avoid E>0 for isochrone for
    // points near the turning points -- perhaps I could change this
    // into an iterative procedure?
    rmin = sqrt(rmin);rmax=sqrt(rmax);
    VecDoub innerg = TargetPot->Forces(in);
    VecDoub outerg = TargetPot->Forces(out);
    double innergrad=0., outergrad=0.;
    for(int j=0;j<3;j++)innergrad-=innerg[j]*in[j];innergrad/=rmin;
    for(int j=0;j<3;j++)outergrad-=outerg[j]*out[j];outergrad/=rmax;
    isoparams_st isos(rmin,rmax,innergrad,outergrad);
    root_find RF(0.01,100);
    double b = RF.findroot(&iso_params,0.,100.,&isos);
    if(b<0.)b=5.;
    double GM = innergrad*(b+sqrt(rmin*rmin+b*b))*(b+sqrt(rmin*rmin+b*b))*sqrt(rmin*rmin+b*b)/rmin;
    // now catch those that give too high an energy
    Isochrone Iso(GM,b);
    // std::vector<int> failures;
    // for(unsigned i=0;i<orbit_samples.size();i++)
    //     if(Iso.H(orbit_samples[i])>0.)
    //         failures.push_back(i);
    // int still_fail=1;
    // while(still_fail){
    //     GM*=1.001;
    //     Iso.set({GM});
    //     still_fail=0;
    //     for(unsigned i=0;i<orbit_samples.size();i++)
    //         if(Iso.H(orbit_samples[i])>0.){ still_fail=1; break;}
    // }
    return {GM,b};
}

VecDoub Actions_Genfunc::actions(const VecDoub &x, void *params){
    // for(auto i: AF->actions(x)) std::cout<<i<<" ";
    // Set the parameters -- N_T, N_max, Total_T
	Actions_Genfunc_data_structure *AA;
    if(params!=nullptr) AA = (Actions_Genfunc_data_structure *)params;
    else{
        double timescale = TargetPot->torb(x);
        AA = new Actions_Genfunc_data_structure(timescale*8., 300, 8, 1e-8,1200,timescale*32.,symmetry=="axisymmetric"?true:false);
    }
    // superlong: 800, 20000, 80000, 3200

	// std::cout<<AA->total_T<<" ";
 //    for(auto i: AF->actions(x)) std::cout<<i<<" ";
    // Integrate the orbit
	Orbit orbit(TargetPot,AA->orbit_eps);
	orbit.integrate(x,AA->total_T,AA->stepsize,false);

    // SymplecticOrbit orbit(TargetPot,50.);
    // orbit.integrate(x,AA->total_T,AA->stepsize,false);

    // Assess the angular momentum of the orbit and choose an appropriate
    // toy potential
	std::vector<int> loopi = loop(orbit.results());
	int total_loop=0;
    for(int i=0;i<3;i++)total_loop+=abs(loopi[i]);

    if(total_loop==0){
		// VecDoub Om = find_box_params_minvar(orbit.results());
        VecDoub Om = find_box_params(orbit.results());
		ToyAct = new Actions_HarmonicOscillator(Om);
	}
	else if(total_loop==1){
        // VecDoub ff = find_isochrone_params_minvar(orbit.results());
        VecDoub ff = find_isochrone_params(orbit.results());
		ToyAct = new Actions_Isochrone(ff[0],ff[1]);
	}
    else if(AA->total_T<AA->maxtimescale){
        AA->total_T*=2.;
        Actions_Genfunc_data_structure AA2 = *AA;
        return actions(x,&AA2);
    }
    else{
        if(debug_genfunc){
            std::cerr<<"Cannot find ascertain orbit type for: ";
            for(auto i:x) std::cerr<<i<<" ";
            std::cerr<<". Assuming box.\n";
        }
        VecDoub Om = find_box_params(orbit.results());
        ToyAct = new Actions_HarmonicOscillator(Om);
        // VecDoub ff = find_isochrone_params(orbit.results());
        // ToyAct = new Actions_Isochrone(ff[0],ff[1]);
    }

	// Fill an array with the appropriate vectors n
    std::vector<std::vector<int>> n_vectors;
    int symNx = 2; if(loopi[0]!=0 or loopi[2]!=0) symNx = 1;
    for(int i=-AA->N_matrix;i<AA->N_matrix+1;i+=symNx)
    // specialise for axisymmetric
    for(int j=(AA->axisymmetric?0:-AA->N_matrix);j<(AA->axisymmetric?0:AA->N_matrix)+1;j+=2)
    for(int k=0;k<AA->N_matrix+1;k+=2){
    	if(i==0 and j==0 and k==0) continue;
   		if (k>0 or (k==0 and j>0) or (k==0 and j==0 and i>0))
   		if(sqrt((double)(i*i+j*j+k*k))<=AA->N_matrix){
    			n_vectors.push_back({i,j,k});
    	}
	}

    // Now initialise the arrays. As a is symmetric we use the packed storage
    // format. We use upper-triangular matrix (uplo='U') below so
    // a11 a12 a13 a14
    //     a22 a23 a24
    //         a33 a34     (aij = aji)
    //             a44

    int N_SIZE = n_vectors.size()+3;
    VecDoub a(N_SIZE*(N_SIZE+1)/2,0),b(N_SIZE,0.);

    int orb_SIZE = orbit.results().size();
    // Upper 3 by 3 of a is identity matrix
    a[0]=orb_SIZE;
    a[2]=orb_SIZE;
    a[5]=orb_SIZE;

    // We iterate over the orbit integration points, calculating the toy
    // actions and adding the appropriate terms to the grid
    VecDoub acts(3,0), angs(3,0);
    VecDoub av_acts(3,0);
    VecDoub cd(N_SIZE-3,0.);
    bool longer_int_window = false, finer_sampling = false;
    VecDoub minAdot(N_SIZE-3,1e7);
    VecDoub maxAdot(N_SIZE-3,-1e7);
    VecDoub angs_pp(3,0.),angs_p(3,0.),angs_u=angs_p, acts_p(3,0.),angs_0(3,0.),angs_pu(3,0.);
    VecDoub dangs(3,0.);

    unsigned nn=0; VecDoub sign = {1.,1.,1.};
    for(auto i: orbit.results()){
        orient_orbit(i,loopi);
        if(total_loop==1)
            sign[1] = SIGN(TargetPot->Lz(i));

    	acts = ToyAct->actions(i);
    	angs = ToyAct->angles(i);
        if(nn>0)
            for(int j=0;j<3;++j)
                if((angs[j]-angs_p[j]+PI/2.*sign[j])*sign[j]<0.)
                    angs_u[j]+=sign[j]*2.*PI;
        angs = angs+angs_u;
        // av_acts = av_acts+acts;
        if(nn==0) angs_0=angs;
        if(nn==1)
            dangs = angs-angs_pu;
        else
            dangs = angs-angs_pp;
        dangs = dangs*0.5;
        if(nn>=1){
            for(int i=0;i<3;i++)dangs[i] = acts[i]*dangs[i];
            av_acts = av_acts+dangs;
        }
        if(nn==AA->N_T-1){
            dangs = angs-angs_pu;
            dangs = dangs*0.5;
            for(int i=0;i<3;i++)dangs[i] = acts[i]*dangs[i];
            av_acts = av_acts+dangs;
        }
        // b = (J,2(n.\theta)cos(n.\theta))
    	for(int j=0;j<3;j++) b[j]+=acts[j];
    	for(int j=0;j<N_SIZE-3;j++){
            double dd = dot_product_int<double>(n_vectors[j],angs);
            cd[j]=cos(dd);
        	b[j+3]+=2.*dot_product_int<double>(n_vectors[j],acts)*cd[j];
            if(dd<minAdot[j])minAdot[j]=dd;
            if(dd>maxAdot[j])maxAdot[j]=dd;
        }
        // Each column after 3rd of a (indexed with j=0,1,...)
        // top three elements = 2 n_j cos(n.\theta) i.e. jth n vector
        // all others (k=0,1,...) = 4 dot(n_j,n_k) cos(n_j.\theta) cos(n_k.\theta) i.e. jth and kth n_vector
        for(int j=0;j<N_SIZE-3;j++){
    		for(int k=0;k<3;k++) a[6+k+j*(j+7)/2]+=2.*n_vectors[j][k]*cd[j];
    		for(int k=0;k<j+1;k++){
    			a[9+k+j*(j+7)/2]+=4.*dot_product<int>(n_vectors[j],n_vectors[k])*cd[j]*cd[k];
    		}
    	}
        acts_p=acts;
        angs_pp = angs_pu;
        angs_pu = angs;
        angs_p = angs-angs_u;
        nn++;
    }

    for(int i=0;i<3;++i)av_acts[i]/=(angs[i]-angs_0[i]);

    double mmA;
    for(int i=0;i<N_SIZE-3;++i){
        mmA =fabs(maxAdot[i]-minAdot[i]);
        if(mmA<2.*PI)
            longer_int_window=true;
        if(mmA/(double)AA->N_T>PI)
            finer_sampling=true;
        if(finer_sampling and longer_int_window) break;
    }

    // Now we use LAPACK to solve
    char uplo = 'U'; // shape of matrix a ==> Upper triangular
    int nrhs = 1;    // width of b
    int LDB = N_SIZE;// length of b
    int info;        // error report -- info=0 ==> successful
    int ipiv[N_SIZE];// indicates what has been switched around in factorize
    dspsv_(&uplo, &N_SIZE, &nrhs, & *a.begin(), ipiv, & *b.begin(), &LDB,&info);

    // Now check the results

    // Add a little check to see if actions differ from their average by more than the greatest average action
    bool fac_of_two = false;
    // for(int i=0;i<3;i++) av_acts[i]/=orbit.results().size();
    // for(int i=0;i<3;i++)
    //     if(fabs(b[i]-av_acts[i])>Max<double>(av_acts)) fac_of_two=true;

    // sort out loop orientations
    if(loopi[0]!=0 and total_loop==1){
        double tmp = av_acts[1];av_acts[1]=av_acts[2];av_acts[2]=tmp;
        tmp = b[1];b[1]=b[2];b[2]=tmp;
    }
    if(loopi[0]*b[2]>0.){ b[2]*=-1.; av_acts[1]*=-1.;}
    if(loopi[2]*b[1]<0.){ b[1]*=-1.; av_acts[2]*=-1.;}


    bool need_to_repeat=false, fatal=false;
    if(longer_int_window){
        AA->total_T*=2.; need_to_repeat=true;
    }
    if(finer_sampling){
        AA->N_T*=2.; need_to_repeat=true;
    }
    if(
      // Check if box and all positive
      (total_loop==0 and std::any_of(b.begin(),b.begin()+3,[](double i){return i<0;}))
        // Check if radial action positive
        or b[0]<0.
        // Check estimates are roughly right
        or fac_of_two
        ){
            AA->total_T*=2.;
            AA->N_T*=2.;
            fatal=true;
            need_to_repeat=true;
    }

    if(debug_genfunc){
        std::cerr<<"Finer sampling: "<<finer_sampling<<
                   ", Longer window: "<<longer_int_window<<
                   ", fatal: "<<fatal<<std::endl;
    }

    if(fatal)
        if(AA->N_T>=AA->NTMAX and AA->total_T>=AA->maxtimescale){
            if(debug_genfunc){
                std::cerr<<"Returning average actions\n";
                std::cerr<<"Returning with "<<AA->N_T/AA->NTMAX
                 <<" fraction of max steps and"<<AA->total_T/AA->maxtimescale
                 <<" fraction of max time"<<std::endl;
            }
            if(symmetry=="triaxial"){
            if(std::any_of(loopi.begin(),loopi.end(),[](int i){return i!=0;}))
                av_acts[0]*=2.;
            if(loopi[0]*loopi[0]==1 and la_switch==true){
                double r2 = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
                double acts;
                if(AF)
                    acts = AF->actions(x)[3];
                else{
                    Actions_TriaxialStackel_Fudge ActS(TargetPot,-r2/10.,-r2/20.);
                    acts = ActS.actions(x)[3];
                }
                if(acts==2){
                    double tmp = av_acts[0];av_acts[0]=av_acts[1];av_acts[1]=tmp;
                }
            }
            }
            // for(auto i: AF->actions(x)) std::cout<<i<<" ";
            return av_acts;
        }
    if(need_to_repeat)
        if(AA->N_T<AA->NTMAX and AA->total_T<AA->maxtimescale){
            Actions_Genfunc_data_structure AA2 = *AA;
            return actions(x,&AA2);
        }

    // If we are here the results are OK or we have reached our limit with no
    // big errors

    // scale radial action
    if(symmetry=="triaxial"){
        if(std::any_of(loopi.begin(),loopi.end(),[](int i){return i!=0;}))
            b[0]*=2.;
    // Now check if inner or outer long-axis loop
    // Is there a clever way of doing this?
    // Let's use the Stackel fudge
    if(loopi[0]*loopi[0]==1 and la_switch==true){
        double r2 = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
        double acts;
        if(AF)
            acts = AF->actions(x)[3];
        else{
            Actions_TriaxialStackel_Fudge ActS(TargetPot,-r2/10.,-r2/20.);
            acts = ActS.actions(x)[3];
        }
        if(acts==2){
            double tmp = b[0];b[0]=b[1];b[1]=tmp;
        }
    }
    }
    if(debug_genfunc)
        std::cerr<<"Returning with "<<AA->N_T/AA->NTMAX
                 <<" fraction of max steps and "<<AA->total_T/AA->maxtimescale
                 <<" fraction of max time"<<std::endl;
    // for(auto i: AF->actions(x)) std::cout<<i<<" ";
	return {b[0],b[1],b[2]};
}

VecDoub Actions_Genfunc::angles(const VecDoub &x, void *params){

    // Set the parameters -- N_T, N_max, Total_T
    Actions_Genfunc_data_structure *AA;
    if(params!=nullptr) AA = (Actions_Genfunc_data_structure *)params;
    else{
        double timescale = TargetPot->torb(x);
        AA = new Actions_Genfunc_data_structure(timescale*8., 300, 8, 1e-8,1200,timescale*32.,symmetry=="axisymmetric"?true:false);
    }
    // Integrate the orbit
    Orbit orbit(TargetPot,AA->orbit_eps);
    orbit.integrate(x,AA->total_T,AA->stepsize,false);

    // SymplecticOrbit orbit(TargetPot,50.);
    // orbit.integrate(x,AA->total_T,AA->stepsize,false);

    // Assess the angular momentum of the orbit and choose an appropriate
    // toy potential
    std::vector<int> loopi = loop(orbit.results());
    int total_loop=0;
    for(int i=0;i<3;i++)total_loop+=abs(loopi[i]);

    if(total_loop==0){
        // VecDoub Om = find_box_params_minvar(orbit.results());
        VecDoub Om = find_box_params(orbit.results());
        ToyAct = new Actions_HarmonicOscillator(Om);
    }
    else if(total_loop==1){
        // VecDoub ff = find_isochrone_params_minvar(orbit.results());
        VecDoub ff = find_isochrone_params(orbit.results());
        ToyAct = new Actions_Isochrone(ff[0],ff[1]);
    }
    else if(AA->total_T<AA->maxtimescale){
        AA->total_T*=2.;
        Actions_Genfunc_data_structure AA2 = *AA;
        return angles(x,&AA2);
    }
    else{
        if(debug_genfunc){
            std::cerr<<"Cannot find ascertain orbit type for: ";
            for(auto i:x) std::cerr<<i<<" ";
            std::cerr<<". Assuming box.\n";
        }
        VecDoub Om = find_box_params(orbit.results());
        ToyAct = new Actions_HarmonicOscillator(Om);
    }

    // Fill an array with the appropriate vectors n
    std::vector<std::vector<int>> n_vectors;
    int symNx = 2; if(loopi[0]!=0 or loopi[2]!=0) symNx = 1;
    for(int i=-AA->N_matrix;i<AA->N_matrix+1;i+=symNx)
    // specialise for axisymmetric
    for(int j=(AA->axisymmetric?0:-AA->N_matrix);j<(AA->axisymmetric?0:AA->N_matrix)+1;j+=2)
    for(int k=0;k<AA->N_matrix+1;k+=2){
        if(i==0 and j==0 and k==0) continue;
        if (k>0 or (k==0 and j>0) or (k==0 and j==0 and i>0))
        if(sqrt((double)(i*i+j*j+k*k))<=AA->N_matrix){
                n_vectors.push_back({i,j,k});
        }
    }

    // Now initialise the arrays. As a is symmetric we use the packed storage
    // format. We use upper-triangular matrix (uplo='U') below so
    // a11 a12 a13 a14
    //     a22 a23 a24
    //         a33 a34     (aij = aji)
    //             a44

    int N_SIZE = 3.*n_vectors.size()+6;
    VecDoub a(N_SIZE*(N_SIZE+1)/2,0),b(N_SIZE,0.);

    int orb_SIZE = orbit.results().size();
    // Upper 3 by 3 of a is identity matrix
    a[0]=orb_SIZE;
    a[2]=orb_SIZE;
    a[5]=orb_SIZE;
    double sum_timeseries = (orb_SIZE-1)*orb_SIZE*0.5*AA->stepsize;
    double sum_timeseries2 = (orb_SIZE-.5)*(orb_SIZE-1)*orb_SIZE/3.*AA->stepsize*AA->stepsize;
    a[6]=sum_timeseries;
    a[11]=sum_timeseries;
    a[17]=sum_timeseries;
    a[9]=sum_timeseries2;
    a[14]=sum_timeseries2;
    a[20]=sum_timeseries2;

    // We iterate over the orbit integration points, calculating the toy
    // actions and adding the appropriate terms to the grid
    VecDoub acts(3,0), angs(3,0);
    VecDoub av_acts(3,0);
    VecDoub sd(n_vectors.size(),0.);
    bool longer_int_window = false, finer_sampling = false;
    VecDoub minAdot(n_vectors.size(),1e7);
    VecDoub maxAdot(n_vectors.size(),-1e7);
    VecDoub angs_pp(3,0.),angs_p(3,0.),angs_u=angs_p, acts_p(3,0.),angs_0(3,0.),angs_pu(3,0.);
    double tm;
    unsigned nn=0; VecDoub sign = {1.,1.,1.};
    for(auto i: orbit.results()){
        tm=AA->stepsize*(double)nn;
        orient_orbit(i,loopi);
        if(total_loop==1)
            sign[1] = SIGN(TargetPot->Lz(i));
        angs = ToyAct->angles(i);
        if(nn>0)
            for(int j=0;j<3;++j)
                if((angs[j]-angs_p[j]+PI/2.*sign[j])*sign[j]<0.)
                    angs_u[j]+=sign[j]*2.*PI;
        angs = angs+angs_u;
        // b = (J,2(n.\theta)cos(n.\theta))
        for(int j=0;j<3;j++) b[j]+=angs[j];
        for(int j=3;j<6;j++) b[j]+=angs[j-3]*tm;
        for(unsigned j=0;j<n_vectors.size();j++){
            double dd = dot_product_int<double>(n_vectors[j],angs);
            sd[j]=2.*sin(dd);
            b[j+6]-=angs[0]*sd[j];
            b[j+6+n_vectors.size()]-=angs[1]*sd[j];
            b[j+6+2.*n_vectors.size()]-=angs[2]*sd[j];
            if(dd<minAdot[j])minAdot[j]=dd;
            if(dd>maxAdot[j])maxAdot[j]=dd;
        }
        // Each column after 3rd of a (indexed with j=0,1,...)
        // top three elements = 2 n_j cos(n.\theta) i.e. jth n vector
        // all others (k=0,1,...) = 4 dot(n_j,n_k) cos(n_j.\theta) cos(n_k.\theta) i.e. jth and kth n_vector
        for(unsigned j=0;j<n_vectors.size();j++){
            a[(j+6)*(j+7)/2]-=sd[j];
            a[1+(j+n_vectors.size()+6)*(j+n_vectors.size()+7)/2]-=sd[j];
            a[2+(j+2.*n_vectors.size()+6)*(j+2.*n_vectors.size()+7)/2]-=sd[j];
            a[3+(j+6)*(j+7)/2]-=tm*sd[j];
            a[4+(j+n_vectors.size()+6)*(j+n_vectors.size()+7)/2]-=tm*sd[j];
            a[5+(j+2.*n_vectors.size()+6)*(j+2.*n_vectors.size()+7)/2]-=tm*sd[j];

            for(int k=0;k<j+1;k++){
                a[6+k+(j+6)*(j+7)/2]+=sd[j]*sd[k];
                a[6+n_vectors.size()+k+(j+n_vectors.size()+6)*(j+n_vectors.size()+7)/2]+=sd[j]*sd[k];
                a[6+2*n_vectors.size()+k+(j+2*n_vectors.size()+6)*(j+2*n_vectors.size()+7)/2]+=sd[j]*sd[k];
            }
        }
        // copy angles and freqs
        acts_p=acts;
        angs_pp = angs_pu;
        angs_pu = angs;
        angs_p = angs-angs_u;
        nn++;
    }

    double mmA;
    for(unsigned i=0;i<n_vectors.size();++i){
        mmA =fabs(maxAdot[i]-minAdot[i]);
        if(mmA<2.*PI)
            longer_int_window=true;
        if(mmA/(double)AA->N_T>PI)
            finer_sampling=true;
        if(finer_sampling and longer_int_window) break;
    }

    // Now we use LAPACK to solve
    char uplo = 'U'; // shape of matrix a ==> Upper triangular
    int nrhs = 1;    // width of b
    int LDB = N_SIZE;// length of b
    int info;        // error report -- info=0 ==> successful
    int ipiv[N_SIZE];// indicates what has been switched around in factorize
    dspsv_(&uplo, &N_SIZE, &nrhs, & *a.begin(), ipiv, & *b.begin(), &LDB,&info);

    // Now check the results

    // sort out loop orientations
    if(loopi[0]!=0 and total_loop==1){
        double tmp = b[1];b[1]=b[2];b[2]=tmp;
        tmp = b[4];b[4]=b[5];b[5]=tmp;
    }
    // if(loopi[0]*b[2]>0.) b[2]*=-1.;
    // if(loopi[2]*b[1]<0.) b[1]*=-1.;


    bool need_to_repeat=false, fatal=false;
    if(longer_int_window){
        AA->total_T*=2.; need_to_repeat=true;
    }
    if(finer_sampling){
        AA->N_T*=2.; need_to_repeat=true;
    }

    if(debug_genfunc){
        std::cerr<<"Finer sampling: "<<finer_sampling<<
                   ", Longer window: "<<longer_int_window<<
                   ", fatal: "<<fatal<<std::endl;
    }

    if(fatal)
        if(AA->N_T>=AA->NTMAX and AA->total_T>=AA->maxtimescale){
            if(debug_genfunc){
                std::cerr<<"Returning average actions\n";
                std::cerr<<"Returning with "<<AA->N_T/AA->NTMAX
                 <<" fraction of max steps and"<<AA->total_T/AA->maxtimescale
                 <<" fraction of max time"<<std::endl;
            }
            return {-99.,-99.,-99.};
        }
    if(need_to_repeat)
        if(AA->N_T<AA->NTMAX and AA->total_T<AA->maxtimescale){
            Actions_Genfunc_data_structure AA2 = *AA;
            return angles(x,&AA2);
        }

    if(debug_genfunc)
        std::cerr<<"Returning with "<<AA->N_T/AA->NTMAX
                 <<" fraction of max steps and "<<AA->total_T/AA->maxtimescale
                 <<" fraction of max time"<<std::endl;
    for(int i=0;i<3;++i){ b[i]=fmod(b[i],TPI); if(b[i]<0.) b[i]+=TPI;}
    return {b[0],b[1],b[2],b[3],b[4],b[5]};
}

VecDoub Actions_Genfunc::full_actions(const VecDoub &x,int NT,int Nsamp,int Nmax){
    Actions_Genfunc_data_structure *AA;
    double timescale = TargetPot->torb(x);
    AA = new Actions_Genfunc_data_structure(timescale*NT, Nsamp, Nmax, 1e-8,4*Nsamp,timescale*4.*NT,symmetry=="axisymmetric"?true:false);
    return actions(x,AA);
}

VecDoub Actions_Genfunc::full_angles(const VecDoub &x,int NT,int Nsamp,int Nmax){
    Actions_Genfunc_data_structure *AA;
    double timescale = TargetPot->torb(x);
    AA = new Actions_Genfunc_data_structure(timescale*NT, Nsamp, Nmax, 1e-8,4*Nsamp,timescale*4.*NT,symmetry=="axisymmetric"?true:false);
    return angles(x,AA);
}

VecDoub Actions_Genfunc_Average::actions(const VecDoub &x, void *params){

    // Set the parameters -- N_T, N_max, Total_T
    Actions_Genfunc_data_structure *AA;
    if(params!=nullptr) AA = (Actions_Genfunc_data_structure *)params;
    else{
        double timescale = TargetPot->torb(x);
        AA = new Actions_Genfunc_data_structure(timescale*8., 200, 6, 1e-3,500,32.*timescale,symmetry=="axisymmetric"?true:false);
    }
    // Integrate the orbit
    Orbit orbit(TargetPot,AA->orbit_eps);
    orbit.integrate(x,AA->total_T,AA->stepsize,false);

    // Assess the angular momentum of the orbit and choose an appropriate
    // toy potential
    std::vector<int> loopi = loop(orbit.results());
    int total_loop=0;
    for(int i=0;i<3;i++)total_loop+=abs(loopi[i]);

    if(total_loop==0){
        VecDoub Om = find_box_params(orbit.results());
        ToyAct = new Actions_HarmonicOscillator(Om);
    }
    else if(total_loop==1){
        VecDoub ff = find_isochrone_params(orbit.results());
        ToyAct = new Actions_Isochrone(ff[0],ff[1]);
    }
    else if(AA->N_T<200000){
        AA->total_T*=10.;
        AA->N_T*=10.;
        Actions_Genfunc_data_structure AA2 = *AA;
        return actions(x,&AA2);
    }
    else{
        std::cerr<<"Cannot find actions for: ";printVector(x);
        VecDoub Om = find_box_params(orbit.results());
        ToyAct = new Actions_HarmonicOscillator(Om);
    }

    // We iterate over the orbit integration points, calculating the toy
    // actions and adding the appropriate terms to the grid
    VecDoub acts(3,0), angs(3,0);
    VecDoub av_acts(3,0),angs_pp(3,0.),angs_p(3,0.),angs_u=angs_p, acts_p(3,0.),angs_0(3,0.),angs_pu(3,0.),dangs(3,0.),sign(3,1.);

    unsigned k=0;
    for(auto i: orbit.results()){
        orient_orbit(i,loopi);
        if(total_loop==1)
            sign[1] = SIGN(TargetPot->Lz(i));
        acts = ToyAct->actions(i);
        angs = ToyAct->angles(i);

        if(k>0)
            for(int j=0;j<3;++j)
                if((angs[j]-angs_p[j]+PI/2.*sign[j])*sign[j]<0.)
                    angs_u[j]+=sign[j]*2.*PI;
        angs = angs+angs_u;
        // av_acts = av_acts+acts;
        if(k==0) angs_0=angs;
        if(k==1)
            dangs = angs-angs_pu;
        else
            dangs = angs-angs_pp;
        dangs = dangs*0.5;
        if(k>=1){
            for(int i=0;i<3;i++)dangs[i] = acts[i]*dangs[i];
            av_acts = av_acts+dangs;
        }
        if(k==AA->N_T-1){
            dangs = angs-angs_pu;
            dangs = dangs*0.5;
            for(int i=0;i<3;i++)dangs[i] = acts[i]*dangs[i];
            av_acts = av_acts+dangs;
        }
        acts_p=acts;
        angs_pp = angs_pu;
        angs_pu = angs;
        angs_p = angs-angs_u;
        k++;
    }

    for(int i=0;i<3;++i)av_acts[i]/=(angs[i]-angs_0[i]);

    if(loopi[0]!=0 and total_loop==1){
        double tmp = av_acts[1];av_acts[1]=av_acts[2];av_acts[2]=tmp;
    }
    if(loopi[0]*av_acts[2]>0.) av_acts[2]*=-1.;
    if(loopi[2]*av_acts[1]<0.) av_acts[1]*=-1.;

    if(((
       // Check if box and all positive
       std::all_of(loopi.begin(),loopi.end(),[](int i){return i==0;})
       and std::any_of(av_acts.begin(),av_acts.end(),[](int i){return i<0;}))
        // Check if radial action positive
        or av_acts[0]<0.)
        // after a certain point give up
        and AA->N_T<200000){
        AA->total_T*=10.;
        AA->N_T*=10.;
        Actions_Genfunc_data_structure AA2 = *AA;
        return actions(x,&AA2);
    }

    // scale radial action
    if(symmetry=="triaxial"){
        if(std::any_of(loopi.begin(),loopi.end(),[](int i){return i!=0;}))
            av_acts[0]*=2.;
    // Now check if inner or outer long-axis loop
    // Is there a clever way of doing this?
    // Let's use the Stackel fudge
    if(loopi[0]*loopi[0]==1 and la_switch==true){
        double r2 = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
        double acts;
        if(AF)
            acts = AF->actions(x)[3];
        else{
            Actions_TriaxialStackel_Fudge ActS(TargetPot,-r2/10.,-r2/20.);
            acts = ActS.actions(x)[3];
        }
        if(acts==2){
            double tmp = av_acts[0];av_acts[0]=av_acts[1];av_acts[1]=tmp;
        }
    }}
    return av_acts;

}
/*int main(int argc,char*argv[]){
    VecDoub X = {1.,.05,1.,0.2,0.7,0.2};
    // StackelTriaxial Iso(1.,-30.,-20.);
    MultipoleExpansion_Triaxial Iso(argv[1]);
    // Orbit O(&Iso);
    // O.integrate(X,10.*Iso.torb(X),.01*Iso.torb(X));
    // O.plot(0,1);
    // HarmonicOscillator Iso({10.,7.,4.});
    // Isochrone Iso(2e6,6.);
    // Actions_Isochrone AI(2e6,6.);
    // printVector(AI.actions(X));
    Actions_Genfunc ActI(&Iso);
    Actions_TriaxialStackel_Fudge ActS(&Iso,-5.,-2.);
    // for(double x = 0.;x<150.;x+=5.){
    // X[3]=x;
    // std::cout<<x<<" ";
    // for(auto i:ActS.actions(X)) std::cout<<i<<" ";
    // for(auto i:ActI.actions(X)) std::cout<<i<<" ";
    // std::cout<<std::endl;
    // }
    // Orbit O(&Iso,1e-3);
    // VecDoub Q = O.integrate(X,5.,0.01);
    // // printVector(ActI.actions(Q,&ASGD));
    // Actions_TriaxialStackel ActS(&Iso);
    printVector(ActS.actions(X));
    printVector(ActI.actions(X));
    // printVector(ActS.actions(Q,&ASGD));
    return 0;
}
*/

VecDoub PhaseSpacePoint(VecDoub actions, VecDoub angles, Potential_JS *Pot){
	Torus T; Actions Acts; Angles Angs;
	WrapperTorusPotential TPot(Pot);
	vec2torus(actions,Acts);
	vec2torus(angles,Angs);
	T.AutoFit(Acts,&TPot,1e-7,700,300,15,5,24,200,24,0);
	return torusPSPT2cartvec(T.Map3D(Angs));
}
// #include "LogPot.h"
// int main(){
// 	double kk = 977.7775;
// 	Logarithmic Log(220.,1.,0.9);
// 	Actions_Genfunc GG(&Log,"axisymmetric",nullptr);
// 	VecDoub X = {10.,4.,7.,80.,260.,40.};
// 	VecDoub aa = GG.actions(X);
// 	VecDoub ang = GG.angles(X);
// 	Torus T;
// 	LogPot Log2(220./kk,0.9,1e-5,0.);
// 	Actions AA;AA[0]=aa[0];AA[1]=aa[2];AA[2]=aa[1];
// 	AA/=kk;
// 	T.AutoFit(AA,&Log2,1e-7,700,300,15,5,24,200,24,0);
// 	Angles Ang;Ang[0]=ang[0];Ang[2]=ang[1];Ang[1]=ang[2];
// 	PSPT FF = T.Map3D(Ang);
// 	for(int i=3;i<6;++i)FF[i]*=kk;
// }

//      VecDoub X = {1e-5,2.,1e-5,0.001,0.1,0.23};
//      // VecDoub X =  {1e-5,1e-5,4.4,70.01,15.1,15.01};//ILAL
//      // StackelTriaxial Pot(3.61*500.,-30.,-20.);
//      // VecDoub X = {26.,0.1,0.1,0.,122.,83.1};
//      // NFW Pot(1.,1.,0.95,0.85);
//      MultipoleExpansion_Triaxial Pot("/data/jls/self_con/June15/iso_ap1.3_az1.6.me");
//      std::cout<<Pot.H(X)<<std::endl;
//      Orbit O(&Pot);
//      double torb = Pot.torb(X);
//      std::cout<<torb<<std::endl;VecDoub Y = O.integrate(X,10.*torb,0.01*torb);
//      O.plot(0,1);
//      Actions_TriaxialStackel_Fudge ActS(&Pot,-5.,-2.);
//      // lmn_orb ActS(&Pot);
//      Actions_Genfunc ActI(&Pot,"triaxial",nullptr);
//      printVector(ActS.actions(Y));
//      printVector(ActI.actions(Y));
//      std::ofstream outfile;outfile.open("tmp");
//      for(auto i: O.results()){
//         for(auto j:ActS.actions(i))outfile<<j<<" ";
//         outfile<<std::endl;
//     }
//     outfile.close();
//     return 0;
//  }
