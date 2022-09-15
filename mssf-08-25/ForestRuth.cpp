#include <iostream>
#include <cmath>
#include "vector.h"

using namespace std;

//Constantes globales
const double GM=1;

//Constantes Forest-Roth
const double Theta=1/(2-pow(2.0,1.0/3));
const double ThetaU2 = Theta/2;
const double UmThetaU2 = (1 - Theta)/2;
const double Um2Theta = 1-2*Theta;

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
    void Arranque(double dt);
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
    tau=r*normF/r.norm();
}

void Cuerpo::Arranque(double dt){
    omega-=tau*(dt/(2*m));
}

void Cuerpo::Mueva_theta(double dt, double coef){
    r+= omega*(dt*coef);
}
void Cuerpo::Mueva_omega(double dt, double coef){
    omega+= tau*(dt*coef/m);
}

int main(){
    Cuerpo Pendulo;
    double t, dt=0.1;
    double omega,  T;
    double r0=5, v0, m0=1;

    omega=sqrt(GM/(r0*r0*r0));
    v0= omega*r0;
    T=2*M_PI/omega;

    //Debug
    //cout<<T<<" "<<v0<<" "<<omega<<endl;
    //cout<<" "<<Theta<<" "<<ThetaU2<<" "<<Um2Theta<<" "<<UmThetaU2<<endl;

    //----------( x0, y0, omegax0, omegay0, m0, R0);
    Pendulo.Inicie( r0, 0, 0, 0, v0, 0 , m0, 0.15);



    for(t=0; t<T; t+=dt){
        cout<<Pendulo.Getx()<<" "<<Pendulo.Gety()<<endl;
        // Mover por Forest-Roth
        Pendulo.Mueva_theta(dt, ThetaU2);
        Pendulo.CalculeTorque();

        Pendulo.Mueva_omega(dt, Theta);

        Pendulo.Mueva_theta(dt,UmThetaU2);
        Pendulo.CalculeTorque();

        Pendulo.Mueva_omega(dt, Um2Theta);

        Pendulo.Mueva_theta(dt,UmThetaU2);
        Pendulo.CalculeTorque();

        Pendulo.Mueva_omega(dt, Theta);

        Pendulo.Mueva_theta(dt, ThetaU2);
    }
    return 0;
}