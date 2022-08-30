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
    vector3D r,rold, V, F;
    double m,R;
public:
    void Inicie(double x0,double y0,double z0,
                    double Vx0,double Vy0, double Vz0,
                    double m0,double R0);
    void CalculeFuerza(void);
    void Arranque(double dt);
    void Mueva_r(double dt, double coef);
    void Mueva_V(double dt, double coef);
    double Getx(void){return r.x();}
    double Gety(void){return r.y();}
    double Getz(void){return r.z();}
};

void Cuerpo::Inicie(double x0,double y0,double z0,
                    double Vx0,double Vy0, double Vz0,
                    double m0,double R0){
    r.load(x0,y0,z0); V.load(Vx0,Vy0,Vz0); m=m0; R=R0;
}

void Cuerpo::CalculeFuerza(void){
    double normF = -GM * m/r.norm2();
    F=r*normF/r.norm();
}

void Cuerpo::Arranque(double dt){
    V-=F*(dt/(2*m));
}

void Cuerpo::Mueva_r(double dt, double coef){
    r+= V*(dt*coef);
}
void Cuerpo::Mueva_V(double dt, double coef){
    V+= F*(dt*coef/m);
}

int main(){
    Cuerpo Planeta;
    double t, dt=0.1;
    double omega,  T;
    double r0=5, v0, m0=1;

    omega=sqrt(GM/(r0*r0*r0));
    v0= omega*r0;
    T=2*M_PI/omega;

    //Debug
    //cout<<T<<" "<<v0<<" "<<omega<<endl;
    //cout<<" "<<Theta<<" "<<ThetaU2<<" "<<Um2Theta<<" "<<UmThetaU2<<endl;

    //----------( x0, y0, Vx0, Vy0, m0, R0);
    Planeta.Inicie( r0, 0, 0, 0, v0, 0 , m0, 0.15);



    for(t=0; t<T; t+=dt){
        cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<endl;
        // Mover por Forest-Roth
        Planeta.Mueva_r(dt, ThetaU2);
        Planeta.CalculeFuerza();

        Planeta.Mueva_V(dt, Theta);

        Planeta.Mueva_r(dt,UmThetaU2);
        Planeta.CalculeFuerza();

        Planeta.Mueva_V(dt, Um2Theta);

        Planeta.Mueva_r(dt,UmThetaU2);
        Planeta.CalculeFuerza();

        Planeta.Mueva_V(dt, Theta);

        Planeta.Mueva_r(dt, ThetaU2);
    }
    return 0;
}