// ============================================================================
/// \file src/Multipole.cpp
// ============================================================================
/// \author Jason Sanders
/// \date 2014-2015
/// Institute of Astronomy, University of Cambridge (and University of Oxford)
// ============================================================================

// ============================================================================
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// ============================================================================
/// \brief Multipole expansion for general densities
///
/// Main class is MultipoleExpansion that implements the multipole expansion.
/// There are inherited classes that are specific for spherical, axisymmetric
/// , triaxial and up-down symmetric density profiles.
/// Also some test density profiles.
//============================================================================

// #include <Python.h>
#include "coordtransforms.h"
#include "potential.h"
#include "GSLInterface/GSLInterface.h"
#include "Multipole.h"
#include "gnuplot/gnuplot_i.h"
// =====
// NOTE: cuba and openmp don't like each other so must set CUBACORES=0
// =====
// ============================================================================
// Triaxial potential code -- see Section 2.5.3 of BT08
// ============================================================================

static double psi_m_integrand(double m, void *p){
    TriaxialDensity *TD = (TriaxialDensity*)p;
    double a22 = TD->get_a(2)*TD->get_a(2);
    double m2 = a22*(1/m-1.);
    return TD->density_m2(m2)*a22/m/m;
}

static double phi_integrand(double s, void *p){
    phi_integrand_st *P = (phi_integrand_st*)p;
    double tau = 1./s; tau*=tau;
    tau-=1.;
    tau*=P->TP->get_a(2);
    tau*=P->TP->get_a(2);
    double a1 = P->TP->get_a(0)*P->TP->get_a(0)+tau,
           a2 = P->TP->get_a(1)*P->TP->get_a(1)+tau,
           a3 = P->TP->get_a(2)*P->TP->get_a(2)+tau;
    double m = P->x[0]*P->x[0]/a1;
    m += P->x[1]*P->x[1]/a2;
    m += P->x[2]*P->x[2]/a3;
    m*=P->TP->get_a(0)*P->TP->get_a(0);
    return (P->TP->get_psi_inf()-P->TP->psi_m(sqrt(m)))*a3/sqrt(a1*a2);
}

static double force_integrand(double s, void *p){
    phi_integrand_st *P = (phi_integrand_st*)p;
    double tau = 1./s; tau*=tau;
    tau-=1.;
    tau*=P->TP->get_a(2);
    tau*=P->TP->get_a(2);
    double deriv = 2.*P->TP->get_a(0)*P->TP->get_a(0);
    deriv*=P->x[P->index]/(P->TP->get_a(P->index)*P->TP->get_a(P->index)+tau);
    double a1 = P->TP->get_a(0)*P->TP->get_a(0)+tau,
           a2 = P->TP->get_a(1)*P->TP->get_a(1)+tau,
           a3 = P->TP->get_a(2)*P->TP->get_a(2)+tau;
    double m = P->x[0]*P->x[0]/a1;
    m += P->x[1]*P->x[1]/a2;
    m += P->x[2]*P->x[2]/a3;
    m*=P->TP->get_a(0)*P->TP->get_a(0);
    return -P->TP->density_m2(m)*a3*deriv/sqrt(a1*a2);
}


double TriaxialPotential::psi_m(double m){
    // Here we calculate \int_0^m^2 dm^2 rho(m^2)
    // we use a change of variables n = a_2^2/(m^2+a_2^2)
    double a22 = TD->get_a(2)*TD->get_a(2);
    double upper = a22/(m*m+a22);
    // integrator Int(10000);
    // return Int.integrate(&psi_m_integrand,upper,1.,1e-2,TD);
    return GaussLegendreQuad(&psi_m_integrand,upper,1.,TD);
}
double TriaxialPotential::Phi(const VecDoub& x){
    // Here we calculate \int_0^inf dtau f(\tau) (see Eq 2.140 BT08)
    // we use a change of variables n = a_2/sqrt(tau+a_2^2)
    phi_integrand_st P(this, x, -1);
    return -2.*PI*conv::G*get_a(1)/get_a(0)*integrate(&phi_integrand,0.,1.,5e-4,&P);
}

VecDoub TriaxialPotential::Forces(const VecDoub& x){
    VecDoub F = {0.,0.,0.};
    for(int i=0;i<3;i++){
        phi_integrand_st P(this, x, i);
        F[i] = 2.*PI*conv::G*get_a(1)/get_a(0)*GaussLegendreQuad(&force_integrand,0.,1.,&P);
    }
    return F;
}

// ============================================================================
// Multipole expansion
// ============================================================================

MultipoleExpansion::MultipoleExpansion(Density *rho,int NR, int NA, int LMAX1, int MMAX1, double a0, double rmin, double rmax,
    bool axisymmetric, bool triaxial, bool flip, double err,bool vocal)
    :rho(rho), NR(NR), NA(NA), a0(a0),rmin(rmin),rmax(rmax),
     axisymmetric(axisymmetric), triaxial(triaxial), flip(flip), err(err), vocal(vocal){

    NA_phi = axisymmetric?1:NA;
    NA_theta = NA;

    if(LMAX1>0) LMAX = LMAX1;
    else{
        if(flip)LMAX = NA_theta*8.;
        else LMAX = NA_theta*4.;
    }
    if(MMAX1>0) MMAX = MMAX1;
    else MMAX = LMAX;
    if(axisymmetric) MMAX = 0;

    GLtheta = std::unique_ptr<GaussLegendreIntegrator>
            (new GaussLegendreIntegrator(2.*NA_theta-1));
    GLphi = std::unique_ptr<GaussLegendreIntegrator>
            (new GaussLegendreIntegrator(2.*NA_phi-1));

    fill_radial_grid(a0,rmin,rmax);

    rho_grid=std::vector<std::vector<VecDoub>>(2*NA_theta-1,std::vector<VecDoub>(2*NA_phi-1,VecDoub(NR,0.)));
    fill_density_grid();

    Phi_grid=std::vector<std::vector<VecDoub>>(LMAX+1,std::vector<VecDoub>(2.*LMAX+1,VecDoub(NR,0.)));
    dPhi_grid=Phi_grid;
    fillPhigrid();
}

MultipoleExpansion::MultipoleExpansion(const std::string& inFile):err(0.){

    std::ifstream infile; infile.open(inFile);
    if(!infile.is_open())std::cerr<<"Can't open "<<inFile<<std::endl;

    infile >>NR>>NA>>LMAX>>MMAX>>a0>>rmin>>rmax>>axisymmetric>>triaxial>>flip;

    NA_phi = axisymmetric?1:NA;
    NA_theta = NA;

    GLtheta = std::unique_ptr<GaussLegendreIntegrator>
            (new GaussLegendreIntegrator(2.*NA_theta-1));

    GLphi = std::unique_ptr<GaussLegendreIntegrator>
            (new GaussLegendreIntegrator(2.*NA_phi-1));

    fill_radial_grid(a0,rmin,rmax);

    rho_grid=std::vector<std::vector<VecDoub>>(2.*NA_theta-1,std::vector<VecDoub>(2.*NA_phi-1,VecDoub(NR,0.)));
    for(unsigned int i=0;i<rho_grid.size();i++)
        for(unsigned int j=0;j<rho_grid[0].size();j++)
            for(unsigned int k=0;k<rho_grid[0][0].size();k++){
                infile>>rho_grid[i][j][k];
            }

    if(axisymmetric) MMAX = 0;
    Phi_grid=std::vector<std::vector<VecDoub>>(LMAX+1,std::vector<VecDoub>(2.*LMAX+1,VecDoub(NR,0.)));
    dPhi_grid=Phi_grid;
    for(unsigned int i=0;i<Phi_grid.size();i++)
        for(unsigned int j=0;j<Phi_grid[0].size();j++)
            for(unsigned int k=0;k<Phi_grid[0][0].size();k++){
                infile>>Phi_grid[i][j][k];
            }
    for(unsigned int i=0;i<dPhi_grid.size();i++)
        for(unsigned int j=0;j<dPhi_grid[0].size();j++)
            for(unsigned int k=0;k<dPhi_grid[0][0].size();k++){
                infile>>dPhi_grid[i][j][k];
            }
    infile.close();

}


void MultipoleExpansion::fill_radial_grid(double aa0, double rrmin, double rrmax){
    // if loggrid use logarithmic grid else sinh
    a0 = aa0; rmax = rrmax; rmin = rrmin;
    double dmin = asinh(rmin/a0);
    double ra = asinh(rmax/a0)-dmin;
    if(loggrid){
        dmin = log(rmin/a0);
        ra = log(rmax/rmin);
    }
    for(int i=0;i<NR;i++){
        delta_grid.push_back(i*ra/(double)(NR-1)+dmin);
        if(loggrid)
            radial_grid.push_back(a0*exp(delta_grid[i]));
        else
            radial_grid.push_back(a0*sinh(delta_grid[i]));
    }
}

void MultipoleExpansion::add_radial_gridpoint(void){
    delta_grid.push_back(2.*delta_grid[NR-1]-delta_grid[NR-2]);
    if(loggrid)
            radial_grid.push_back(a0*exp(delta_grid.back()));
        else
            radial_grid.push_back(a0*sinh(delta_grid.back()));
    NR++; rmax = radial_grid.back();
    for(unsigned int i=0;i<dPhi_grid.size();i++)
        for(unsigned int j=0;j<dPhi_grid[0].size();j++){
            dPhi_grid[i][j].push_back(0.);Phi_grid[i][j].push_back(0.);
        }
}

void MultipoleExpansion::fill_density_grid(void){

    VecDoub xmin = {-1.,0.}, xmax = {1.,2.*PI};
    if(flip) xmax[0]=0.;
    if(triaxial){ xmax[0]=0.; xmax[1]=PI/2.;}

    double xm=0.5*(xmax[0]+xmin[0]), xr=0.5*(xmax[0]-xmin[0]);
    double ym=0.5*(xmax[1]+xmin[1]), yr=0.5*(xmax[1]-xmin[1]);

    #pragma omp parallel for schedule(dynamic) collapse(2)
    for(int k=0;k<NR;k++){
    for(int i=0;i<NA_theta;i++){
        double dx=xr*GLtheta->abscissa_and_weight_half(i)[0];
        double tu = acos(xm+dx), td = acos(xm-dx);
        for(int j=0;j<NA_phi;j++){
            double dy=yr*GLphi->abscissa_and_weight_half(j)[0];
            double pu = ym-dy, pd = ym+dy;
                rho_grid[NA_theta-1+i][NA_phi-1+j][k]=
                    rho->density_sph({radial_grid[k],pu,tu});
                if(vocal)
                    std::cerr<<OUTPUT(i)<<OUTPUT(j)<<OUTPUT(radial_grid[k])
                <<OUTPUTE(rho_grid[NA_theta-1+i][NA_phi-1+j][k]);
                if(i==0 and j==0) continue;
                rho_grid[NA_theta-1-i][NA_phi-1+j][k]=
                    rho->density_sph({radial_grid[k],pu,td});
                if(j==0) continue;
                rho_grid[NA_theta-1+i][NA_phi-1-j][k]=
                    rho->density_sph({radial_grid[k],pd,tu});
                if(i==0) continue;
                rho_grid[NA_theta-1-i][NA_phi-1-j][k]=
                    rho->density_sph({radial_grid[k],pd,td});
            }
        }
    }
}

void MultipoleExpansion::extend_density_grid(void){
    // add extra radial point to density grid

    VecDoub xmin = {-1.,0.}, xmax = {1.,2.*PI};
    if(flip) xmax[0]=0.;
    if(triaxial){ xmax[0]=0.; xmax[1]=PI/2.;}

    double xm=0.5*(xmax[0]+xmin[0]), xr=0.5*(xmax[0]-xmin[0]);
    double ym=0.5*(xmax[1]+xmin[1]), yr=0.5*(xmax[1]-xmin[1]);

    #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<NA_theta;i++){
        double dx=xr*GLtheta->abscissa_and_weight_half(i)[0];
        double tu = acos(xm+dx), td = acos(xm-dx);
        for(int j=0;j<NA_phi;j++){
            double dy=yr*GLphi->abscissa_and_weight_half(j)[0];
            double pu = ym-dy, pd = ym+dy;
                rho_grid[NA_theta-1+i][NA_phi-1+j].push_back(
                    rho->density_sph({radial_grid.back(),pu,tu}));
                if(vocal)
                    std::cerr<<OUTPUT(i)<<OUTPUT(j)<<OUTPUT(radial_grid.back())
                <<OUTPUTE(rho_grid[NA_theta-1+i][NA_phi-1+j].back());
                if(i==0 and j==0) continue;
                rho_grid[NA_theta-1-i][NA_phi-1+j].push_back(
                    rho->density_sph({radial_grid.back(),pu,td}));
                if(j==0) continue;
                rho_grid[NA_theta-1+i][NA_phi-1-j].push_back(
                    rho->density_sph({radial_grid.back(),pd,tu}));
                if(i==0) continue;
                rho_grid[NA_theta-1-i][NA_phi-1-j].push_back(
                    rho->density_sph({radial_grid.back(),pd,td}));
            }
        }
}

double MultipoleExpansion::rholm_s(int i, int l, int m){
// do integral over cos(theta) and phi

    VecDoub xmin = {-1.,0.}, xmax = {1.,2.*PI};
    double fac = 1.;
    if(flip){ xmax[0]=0.; fac=2.;}
    if(triaxial){ xmax[0]=0.; xmax[1]=PI/2.;
              fac=8.; if(m%4==2) fac*=-1.;}

    double rr=0.0,ss=0.,dx,dy;
    double xm=0.5*(xmax[0]+xmin[0]), xr=0.5*(xmax[0]-xmin[0]);
    double ym=0.5*(xmax[1]+xmin[1]), yr=0.5*(xmax[1]-xmin[1]);
    double tt1, tt2, tt3, tt4;

    for (int j=0; j<NA_theta; j++){
        dx=xr*GLtheta->abscissa_and_weight_half(j)[0];
        ss = 0.;
        for (int k=0; k<NA_phi; k++){
            dy=yr*GLphi->abscissa_and_weight_half(k)[0];
            tt1 = real_spherical_harmonic_2l_1(xm+dx,ym+dy,l,m);
            tt2 = real_spherical_harmonic_2l_1(xm-dx,ym+dy,l,m);
            tt3 = real_spherical_harmonic_2l_1(xm+dx,ym-dy,l,m);
            tt4 = real_spherical_harmonic_2l_1(xm-dx,ym-dy,l,m);
            ss=ss+(k==0?0.5:1.)*GLphi->abscissa_and_weight_half(k)[1]*(
                 tt1*rho_grid[NA_theta-1+j][NA_phi-1+k][i]
                 +tt2*rho_grid[NA_theta-1-j][NA_phi-1+k][i]
                 +tt3*rho_grid[NA_theta-1+j][NA_phi-1-k][i]
                 +tt4*rho_grid[NA_theta-1-j][NA_phi-1-k][i]);
        }
        rr=rr+(j==0?0.5:1.)*GLtheta->abscissa_and_weight_half(j)[1]*ss;
    }
    return fac*rr*xr*yr;
}

void MultipoleExpansion::fillPhigrid(){
    VecDoub rlm(NR,0), y2(NR,0);
    double xl, xh, m_min, m_max, m_skip=1, l_skip=1;

    if(flip or triaxial) l_skip = 2;

    for(int l=0;l<LMAX+1;l+=l_skip){

    m_min = MIN(l,MMAX); m_max = MIN(l,MMAX);
    if(axisymmetric){ m_min=0; m_max=0;}
    if(triaxial){ m_min=0;m_skip = 2;}

    for(int m=-m_min;m<m_max+1;m+=m_skip){

        // 1. Calculate rho*Ylm*r^(l+2) at radial grid points
        #pragma omp parallel for schedule(dynamic)
        for(int i=0;i<NR;i++){
            double drr = radial_grid[i];
            double jac = sqrt(a0*a0+drr*drr);
            if(loggrid)
                jac = drr;
            rlm[i]=rholm_s(i,l,m)
                        *pow(radial_grid[i],l+2)
                        *jac;
        }

        // Make a spline in rlm -- not currently used
        // spline(delta_grid,rlm,1e30,1e30,y2);

        double dx = (delta_grid[1]-delta_grid[0]), P1;//, dx3;
        P1 = 0.;
        // 1.1 First estimate integral from 0 to first grid point
        P1 = 0.5*rlm[0]*radial_grid[0];
        if(loggrid) P1/=radial_grid[0];
        else P1/=sqrt(a0*a0+radial_grid[0]*radial_grid[0]);

        Phi_grid[l][m+l][0]=P1/pow(radial_grid[0],l+1);
        dPhi_grid[l][m+l][0]=-(l+1)*P1/pow(radial_grid[0],l+2);

        // 1.2 Now trapezium rule the rest of the integral
        for(int i=1;i<NR;i++){
            xh = delta_grid[i]; xl=delta_grid[i-1];
            dx = xh-xl; //dx3 = dx*dx*dx;
            P1 +=   0.5*(rlm[i-1]+rlm[i])*dx;
                    //+dx3*(y2[i]+y2[i-1])/12.;
            Phi_grid[l][m+l][i]=P1/pow(radial_grid[i],l+1);
            dPhi_grid[l][m+l][i]=-(l+1)*P1/pow(radial_grid[i],l+2);
        }

        // 2. Calculate rho*Ylm*r^(-l+1) at radial grid points
        for(int i=0;i<NR;i++) rlm[i]=rlm[i]*pow(radial_grid[i],-2*l-1);

        // Make a spline in rlm -- not currently used
        // spline(delta_grid,rlm,1e30,1e30,y2);

        double P2 = 0.;

        // 2.1 Now trapezium rule the rest of the integral
        for(int i=NR-1;i>0;i--){
            xh = delta_grid[i]; xl=delta_grid[i-1];
            dx = xh-xl; //dx3 = dx*dx*dx;
            P2 += 0.5*(rlm[i-1]+rlm[i])*dx;
                    //+dx3*(y2[i]+y2[i-1])/12.;
            Phi_grid[l][m+l][i-1]+=P2*pow(radial_grid[i-1],l);
            dPhi_grid[l][m+l][i-1]+=l*P2*pow(radial_grid[i-1],l-1);
        }

    }}// end of m and l loops

    // Now add grid-points until total mass enclosed is within 1% of
    // target mass in Density class
    if(rho->mass()>0. and err>0.){
        if(fabs(rho->mass()-Mass(radial_grid[NR-1]))>err){
            add_radial_gridpoint();
            extend_density_grid();
            fillPhigrid();
            return;
        }
        if(vocal) std::cout<<"New grid size: NR = "<<NR<<", rmax="<<rmax<<std::endl;
    }

    if(vocal) std::cout<<"Potential calculated: M(r_max)="<<
            	Mass(radial_grid[NR-1])<<std::endl;
    return;
}

double MultipoleExpansion::Phi(const VecDoub& x){
    VecDoub rpt = conv::CartesianToSphericalPolar(x);
    double p = 0., l1; int m_min, m_max, l_skip=1, m_skip=1;

    if(flip or triaxial) l_skip = 2;
    m_min = -MMAX; m_max = MMAX;
    if(axisymmetric){ m_min=0; m_max=0;}
    if(triaxial){ m_min=0;m_skip = 2;}

    for(int m=-m_min;m<m_max+1;m+=m_skip){
        double lm[LMAX-m+1];
        real_spherical_harmonic_2l_1_array(cos(rpt[2]),rpt[1],LMAX,m,lm);
        for(int l=m;l<LMAX+1;l+=l_skip){
            // outside grid quadratic
            if(rpt[0]<radial_grid[0])
                l1 = quad_extrapolate<double>(rpt[0],
                    {radial_grid[0],radial_grid[1],radial_grid[2]},
                    {Phi_grid[l][m+l][0],
                     Phi_grid[l][m+l][1],Phi_grid[l][m+l][2]});
            // outside grid assume vacuum
            else if(rpt[0]>radial_grid[NR-1])
                l1 = Phi_grid[l][m+l][NR-1]*pow(radial_grid[NR-1]/rpt[0],l+1);
            else l1 = linterp(radial_grid,Phi_grid[l][m+l],rpt[0]);
            p+=l1*lm[l-m];
        }
    }
    return -p*conv::FPG;
}

VecDoub MultipoleExpansion::Forces(const VecDoub& x){

    double dPdct = 0., dPdp = 0., dPdr=0.;

    VecDoub rpt = conv::CartesianToSphericalPolar(x);
    double l1; int m_min, m_max, l_skip=1, m_skip=1;

    if(flip or triaxial) l_skip = 2;
    m_min = -MMAX; m_max = MMAX;
    if(axisymmetric){ m_min=0; m_max=0;}
    if(triaxial){ m_min=0;m_skip = 2;}
    if(rpt[0]<radial_grid[0]) m_min=0.;

    for(int m=-m_min;m<m_max+1;m+=m_skip){
        double lm[LMAX-m+1],lmt[LMAX-m+1],lmp[LMAX-m+1];
        real_spherical_harmonic_2l_1_deriv_array(
            cos(rpt[2]),rpt[1],LMAX,m,lm,lmt,lmp);
        for(int l=m;l<LMAX+1;l+=l_skip){
            if(rpt[0]<radial_grid[0])
                l1 = quad_extrapolate<double>(rpt[0],
                    {radial_grid[0],radial_grid[1],radial_grid[2]},
                    {Phi_grid[l][m+l][0],
                     Phi_grid[l][m+l][1],Phi_grid[l][m+l][2]});
            // outside grid assume vacuum
            else if(rpt[0]>radial_grid[NR-1])
                l1 = Phi_grid[l][m+l][NR-1]*pow(radial_grid[NR-1]/rpt[0],l+1);
            else l1 = linterp(radial_grid,Phi_grid[l][m+l],rpt[0]);

            dPdp+=l1*lmp[l-m];
            dPdct+=l1*lmt[l-m];

            // for small r linear
            if(rpt[0]<radial_grid[0]){
                double grad = (dPhi_grid[l][m+l][1]-dPhi_grid[l][m+l][0])/
                           (radial_grid[1]-radial_grid[0]);
                l1 = grad*(rpt[0]-radial_grid[0])+dPhi_grid[l][m+l][0];
            }
            // outside grid assume vacuum
            else if(rpt[0]>radial_grid[NR-1])
                l1 = dPhi_grid[l][m+l][NR-1]*pow(radial_grid[NR-1]/rpt[0],l+2);
            else l1 = linterp(radial_grid,dPhi_grid[l][m+l],rpt[0]);

            dPdr+=-l1*lm[l-m];
            // just use l=0, m=0 for central region
            if(rpt[0]<radial_grid[0])
                break;
            }
        if(rpt[0]<radial_grid[0])
                break;
    }
    VecDoub f = {0.,0.,0.}; double R2 = x[0]*x[0]+x[1]*x[1];
    for(int i=0;i<3;i++) f[i] = dPdr*x[i]/rpt[0];
    f[0]+= x[1]/R2*dPdp+x[0]*x[2]/pow(rpt[0],3.)*dPdct;
    f[1]+= -x[0]/R2*dPdp+x[1]*x[2]/pow(rpt[0],3.)*dPdct;
    f[2]+=-(1.-pow(x[2]/rpt[0],2.))/rpt[0]*dPdct;

    return f*(-conv::FPG);
}

void MultipoleExpansion::output_to_file(const std::string& ofile){
    std::ofstream outfile; outfile.open(ofile);

    outfile <<NR<<" "<<NA<<" "<<LMAX<<" "<<MMAX<<" "<<a0<<" "<<rmin<<" "
            <<rmax<<" "<<axisymmetric<<" "<<triaxial<<" "<<flip<<std::endl;

    for(auto i:rho_grid)
        for(auto j:i){
            for(auto k:j)
                outfile<<k<<" ";
            outfile<<std::endl;
        }
    for(auto i:Phi_grid)
        for(auto j:i){
            for(auto k:j)
                outfile<<k<<" ";
            outfile<<std::endl;
        }
    for(auto i:dPhi_grid)
        for(auto j:i){
            for(auto k:j)
                outfile<<k<<" ";
            outfile<<std::endl;
        }
    outfile.close();
}

void MultipoleExpansion::visualize(const std::string& ofile){
    std::ofstream outfile; outfile.open(ofile);
    VecDoub xmin = {-1.,0.}, xmax = {1.,2.*PI};
    if(flip) xmax[0]=0.;
    if(triaxial){ xmax[0]=0.; xmax[1]=PI/2.;}

    double xm=0.5*(xmax[0]+xmin[0]), xr=0.5*(xmax[0]-xmin[0]);
    double ym=0.5*(xmax[1]+xmin[1]), yr=0.5*(xmax[1]-xmin[1]);

    VecDoub x;
    for(int i=0;i<NA_theta;i++){
        double dx=xr*GLtheta->abscissa_and_weight_half(i)[0];
        double tu = acos(xm+dx), td = acos(xm-dx);
        for(int j=0;j<NA_phi;j++){
            double dy=yr*GLphi->abscissa_and_weight_half(j)[0];
            double pu = ym-dy, pd = ym+dy;
            for(int k=0;k<NR;k++){
                VecDoub r = {radial_grid[k],pu,tu};
                x = conv::SphericalPolarToCartesian(r);
                for(auto i:x)outfile<<i<<" ";for(auto i:r)outfile<<i<<" ";
                outfile<<rho_grid[NA_theta-1+i][NA_phi-1+j][k]<<" "<<Phi(x)<<std::endl;
                if(i==0 and j==0) continue;
                r[2]=td;
                x = conv::SphericalPolarToCartesian(r);
                for(auto i:x)outfile<<i<<" ";for(auto i:r)outfile<<i<<" ";
                outfile<<rho_grid[NA_theta-1-i][NA_phi-1+j][k]<<" "<<Phi(x)<<std::endl;
                if(j==0) continue;
                r[1]=pd;r[2]=tu;
                x = conv::SphericalPolarToCartesian(r);
                for(auto i:x)outfile<<i<<" ";for(auto i:r)outfile<<i<<" ";
                outfile<<rho_grid[NA_theta-1+i][NA_phi-1-j][k]<<" "<<Phi(x)<<std::endl;
                if(i==0) continue;
                r[2]=td;
                x = conv::SphericalPolarToCartesian(r);
                for(auto i:x)outfile<<i<<" ";for(auto i:r)outfile<<i<<" ";
                outfile<<rho_grid[NA_theta-1-i][NA_phi-1-j][k]<<" "<<Phi(x)<<std::endl;
            }
        }
    }
    outfile.close();
}

void test_multipole_spherical(int p=0){
    TestDensity_Hernquist rho(1.,1.,{1.,1.,1.});
    MultipoleExpansion_Spherical ME(&rho,100,1.,0.01,200.);
    VecDoub X = {1e-5,1e-5,1e-5};

    int NMAX = 200;
    VecDoub xx(NMAX,0), exact(NMAX,0), triaxial(NMAX,0), multipole(NMAX,0);

    // #pragma omp parallel for schedule(dynamic)
    for(int xn = 0; xn<NMAX; xn++){
        double x = 0.0001*(double)xn+.0001;
        xx[xn] = x;
        VecDoub X2 = X;
        X2[p]=x;
        exact[xn] = rho.Phi(X2);
        multipole[xn] = ME.Phi(X2);
    }

    Gnuplot G("lines ls 1");
    G.set_xrange(0.9*Min(xx),1.1*Max(xx));
    G.set_yrange(0.9*Min(multipole),1.1*Max(multipole));
    G.set_xlabel("x").set_ylabel("Potential");
    G.savetotex("multipole_density_test").plot_xy(xx,exact);
    G.set_style("lines ls 3").plot_xy(xx,multipole);
    G.outputpdf("multipole_density_test");
}


void test_multipole(int p=0,std::string qq="No"){
    TestDensity_Hernquist rho(1.,10.,{1.,0.9,0.7});
    TriaxialPotential T(&rho,1e6);
    MultipoleExpansion ME2(&rho,100,30,16,6,1.,0.01,2000.,false,true,true);
    ME2.output_to_file("me.tmp");
    ME2.visualize("me.vis");
    MultipoleExpansion ME("me.tmp");
    VecDoub X = {0.001,0.001,0.01};
    double centre2 = 0.;//T.Phi(X);
    double centre3 = 0.;//ME2.Phi(X);

    int NMAX = 200;
    VecDoub xx(NMAX,0), exact(NMAX,0), triaxial(NMAX,0), multipole(NMAX,0);

    #pragma omp parallel for schedule(dynamic)
    for(int xn = 0; xn<NMAX; xn++){
        double x = (double)xn+1.;
        xx[xn] = x;
        VecDoub X2=X;
        X2[p]=x;
        if(qq=="general"){
            X[0]=x/2.;X[1]=x/3.;X[2]=x/2.;
        }
        triaxial[xn] = (T.Phi(X2)-centre2);
        multipole[xn] = (ME2.Phi(X2)-centre3);
    }

    Gnuplot G("lines ls 1");
    G.set_xrange(0.9*Min(xx),1.1*Max(xx));
    G.set_yrange(0.9*Min(multipole),1.1*Max(multipole));
    G.set_xlabel("x").set_ylabel("Potential");
    G.savetotex("multipole_density_test").plot_xy(xx,triaxial);
    G.set_style("lines ls 3").plot_xy(xx,multipole);
    G.outputpdf("multipole_density_test");
}

void test_multipole_axisymmetric(int p=0){
    Miyamoto_NagaiDensity rho(1.,1.,0.7);
    MultipoleExpansion_Axisymmetric ME2(&rho,100,30,16,1.,0.01,100.);
    ME2.output_to_file("me.tmp");
    ME2.visualize("me.vis");
    MultipoleExpansion ME("me.tmp");
    VecDoub X = {0.001,0.001,0.01};
    double centre2 = rho.Phi(X);
    double centre3 = ME2.Phi(X);

    int NMAX = 50;
    VecDoub xx(NMAX,0), exact(NMAX,0), triaxial(NMAX,0), multipole(NMAX,0);

    #pragma omp parallel for schedule(dynamic)
    for(int xn = 0; xn<NMAX; xn++){
        double x = 0.1*(double)xn;
        xx[xn] = x;
        VecDoub X2=X;
        X2[p]=x;
        triaxial[xn] = (rho.Phi(X2)-centre2);
        multipole[xn] = (ME2.Phi(X2)-centre3);
    }

    Gnuplot G("lines ls 1");
    G.set_xrange(0.9*Min(xx),1.1*Max(xx));
    G.set_yrange(0.9*Min(multipole),1.1*Max(multipole));
    G.set_xlabel("x").set_ylabel("Potential");
    G.savetotex("multipole_density_test").plot_xy(xx,triaxial);
    G.set_style("lines ls 3").plot_xy(xx,multipole);
    G.outputpdf("multipole_density_test");
}


void test_multipole_stackel(int p=0,std::string qq="No"){

    TestDensity_Stackel rho(1.,-30.,-10.);
    TriaxialPotential T(&rho,1e8);
    MultipoleExpansion ME(&rho,100,12,8,-1,1.,0.01,100.,false,true,true);
    ME.visualize("me.vis");
    VecDoub X = {1.,1.,1.};
    double centre  = rho.Phi(X);
    double centre2 = T.Phi(X);
    double centre3 = ME.Phi(X);

    int NMAX = 100;
    VecDoub xx(NMAX,0), exact(NMAX,0), triaxial(NMAX,0), multipole(NMAX,0);

    #pragma omp parallel for schedule(dynamic)
    for(int xn = 0; xn<NMAX; xn++){
        double x = (double)xn+.1;
        xx[xn] = x;
        VecDoub X2 = X;
        X2[p]=x;
        if(qq=="general"){
            X[0]=x/2.;X[1]=x/2.;X[2]=x/2.;
        }
        exact[xn] = (rho.Phi(X2)-centre);
        triaxial[xn] = (T.Phi(X2)-centre2);
        multipole[xn] = (ME.Phi(X2)-centre3);
    }

    Gnuplot G("lines ls 1");
    G.set_xrange(0.9*Min(xx),1.1*Max(xx));
    G.set_yrange(0.9*Min(exact),1.1*Max(exact));
    G.set_xlabel("x").set_ylabel("Potential");
    G.savetotex("multipole_density_test").plot_xy(xx,exact);
    // G.set_style("lines ls 2").plot_xy(xx,triaxial);
    G.set_style("lines ls 2").plot_xy(xx,multipole);
    G.outputpdf("multipole_density_test");
}

void test_multipole_forces_stackel(int p=0, int q=0,std::string qq="No"){

    TestDensity_Stackel rho(1.,-30.,-10.);
    TriaxialPotential T(&rho,1e8);
    MultipoleExpansion ME(&rho,100,30,16,-1,10.,0.01,100.,false,true,true);
    ME.visualize("me.vis");

    int NMAX = 50;
    VecDoub xx(NMAX,0), exact(NMAX,0), triaxial(NMAX,0), multipole(NMAX,0);

    #pragma omp parallel for schedule(dynamic)
    for(int xn = 0; xn<NMAX; xn++){
        double x = (double)xn+.1;
        xx[xn] = x;
        VecDoub X(3,1);
        X[p]=x;
        if(qq=="general"){
            X[0]=x/2.;X[1]=x/3.;X[2]=x;
        }
        exact[xn] = rho.Forces(X)[q];
        triaxial[xn] = T.Forces(X)[q];
        multipole[xn] = ME.Forces(X)[q];
    }

    Gnuplot G("lines ls 1");
    G.set_xrange(0.9*Min(xx),1.1*Max(xx));
    G.set_yrange(0.9*Min(exact),1.1*Max(exact));
    G.set_xlabel("x").set_ylabel("Potential");
    G.savetotex("multipole_density_test").plot_xy(xx,exact);
    G.set_style("lines ls 2").plot_xy(xx,triaxial);
    G.set_style("lines ls 3").plot_xy(xx,multipole);
    G.outputpdf("multipole_density_test");
}

void test_multipole_forces(int p=0, int q=0, std::string qq = "No"){
    TestDensity_Hernquist rho(1.,10.,{1.,0.6,0.3});
    // TestDensity_Isochrone rho(1.,10.,{1.,0.9,0.7});
    TriaxialPotential T(&rho,1e6);
    MultipoleExpansion ME(&rho,200,20,8,-1,10.,0.001,5000.,false,true,true);

    int NMAX = 100;
    VecDoub xx(NMAX,0), exact(NMAX,0), multipole(NMAX,0), triaxial(NMAX,0);

    // #pragma omp parallel for schedule(dynamic)
    for(int xn = 0; xn<NMAX; xn++){
        double x = 0.01*(double)xn+1e-3;
        xx[xn] = x;
        VecDoub X(3,1e-4);
        X[p]=x;
        if(qq=="general"){
            X[0]=x/2.;X[1]=x;X[2]=x/4.;
        }
        exact[xn] = rho.Forces(X)[q];
        triaxial[xn] = T.Forces(X)[q];
        multipole[xn] = ME.Forces(X)[q];
    }

    Gnuplot G("lines ls 1");
    G.set_xrange(0.9*Min(xx),1.1*Max(xx));
    G.set_yrange(1.1*Min(exact),1e-6);
    G.set_xlabel("x").set_ylabel("Force");
    G.savetotex("multipole_density_test").plot_xy(xx,multipole);
    G.set_style("lines ls 3").plot_xy(xx,triaxial);
    G.outputpdf("multipole_density_test");
}

void multipole_tests(void){

    // std::cout<<"Simple spherical"<<std::endl;
    // test_multipole_spherical();
    // return;

    // std::cout<<"Stackel potential along x, "<<std::endl;
    // test_multipole_stackel(0);
    // std::cout<<"along y, "<<std::endl;
    // test_multipole_stackel(1);
    // std::cout<<"along z."<<std::endl;
    // test_multipole_stackel(2);

    test_multipole_forces_stackel(1,0,"general");

    // std::cout<<"Stackel potential forces along y in x direction, "<<std::endl;
    // test_multipole_forces_stackel(2,2);
    // std::cout<<"in y direction"<<std::endl;
    // test_multipole_forces_stackel(1,1);
    // std::cout<<"in z direction"<<std::endl;
    // test_multipole_forces_stackel(1,2);

    // std::cout<<"Triaxial Hernquist potential along x, "<<std::endl;
    // test_multipole(0);
    // std::cout<<"along y, "<<std::endl;
    // test_multipole(1);
    // std::cout<<"along z."<<std::endl;

    std::cout<<"Triaxial Hernquist potential forces along x in x direction, "<<std::endl;
    // test_multipole_forces(2,2);
    // std::cout<<"in y direction"<<std::endl;
    // test_multipole_forces(1,1);
    // std::cout<<"in z direction"<<std::endl;
    // test_multipole_forces(2,2);
    // std::cout<<"in radial direction"<<std::endl;
    // test_multipole_forces(2,0,"general");

    // std::cout<<"Miyamoto-Nagai potential along x, "<<std::endl;
    // test_multipole_axisymmetric(0);
    // std::cout<<"along y, "<<std::endl;
    // test_multipole_axisymmetric(1);
    // std::cout<<"along z."<<std::endl;
    // test_multipole_axisymmetric(2);

}

// int main(){
//     multipole_tests();
// }

// int main(int argc, char*argv[]){
//     MultipoleExpansion_Triaxial SCM(argv[1]);
//     std::cout<<SCM.torb({0.1,1.,0.1,0.1,0.1,0.7})<<std::endl;
// }
