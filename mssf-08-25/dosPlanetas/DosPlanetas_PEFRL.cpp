#include <iostream>
#include <cmath>
#include "vector.h"

using namespace std;

//Constantes globales
const double GM=1;

//Constantes PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e00;
const double Chi=-0.6626458266981849e-1;

const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

class Cuerpo;

class Cuerpo{
private:
    vector3D r,rold, omega, tau;
    double m,R;
public:
    void Inicie(double x0,double y0,double z0,
                    double omegax0,double omegay0, double omegaz0,
                    double m0,double R0);
    void CalculeTorque(void);
    void Mueva_theta(double dt, double coef);
    void Mueva_omega(double dt, double coef);
    double Getx(void){return r.x();}
    double Gety(void){return r.y();}
    double Getz(void){return r.z();}
};

void Cuerpo::Inicie(double x0,double y0,double z0,
                    double omegax0,double omegay0, double omegaz0,
                    double m0,double R0){
    r.load(x0,y0,z0); omega.load(omegax0,omegay0,omegaz0); m=m0; R=R0;
}

void Cuerpo::CalculeTorque(void){
    double normF = -GM * m/r.norm2();
    tau=r*(normF/r.norm());
}

void Cuerpo::Mueva_theta(double dt, double coef){
    r+= omega*(dt*coef);
}
void Cuerpo::Mueva_omega(double dt, double coef){
    omega+= tau*(dt*coef/m);
}

int main(){
    Cuerpo Pendulo;
    double t, dt=0.001;
    double r0=5, v0, m0=1;
    double omega, T;
    
    omega=sqrt(GM)*pow(r0,-1.5);
    v0= omega*r0;
    T=2*M_PI/omega;

    //Debug
    //cout<<T<<" "<<v0<<" "<<omega<<endl;
    //cout<<" "<<Zeta<<" "<<" "<<Chi<<" "<<" "<<Lambda<<" "<<Coeficiente1<<" "<<Coeficiente2<<endl;

    //------------( x0,y0,z0,omegax0,omegay0,omegaz0, m0, R0);
    Pendulo.Inicie( r0, 0, 0,  0, 0.75*v0,  0, m0, 0.5);



    for(t=0; t<T; t+=dt){
        cout<<Pendulo.Getx()<<" "<<Pendulo.Gety()<<endl;
        // Mover por PEFRL
        Pendulo.Mueva_theta(dt,Zeta);

        Pendulo.CalculeTorque();
        Pendulo.Mueva_omega(dt, Coeficiente1);

        Pendulo.Mueva_theta(dt, Chi);

        Pendulo.CalculeTorque();
        Pendulo.Mueva_omega(dt, Lambda);

        Pendulo.Mueva_theta(dt, Coeficiente2);

        Pendulo.CalculeTorque();
        Pendulo.Mueva_omega(dt, Lambda);

        Pendulo.Mueva_theta(dt, Chi);

        Pendulo.CalculeTorque();
        Pendulo.Mueva_omega(dt, Coeficiente1);

        Pendulo.Mueva_theta(dt,Zeta);
    }
    return 0;
}