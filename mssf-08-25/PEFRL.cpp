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
    vector3D r,rold, V, F;
    double m,R;
public:
    void Inicie(double x0,double y0,double z0,
                    double Vx0,double Vy0, double Vz0,
                    double m0,double R0);
    void CalculeFuerza(void);
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
    F=r*(normF/r.norm());
}

void Cuerpo::Mueva_r(double dt, double coef){
    r+= V*(dt*coef);
}
void Cuerpo::Mueva_V(double dt, double coef){
    V+= F*(dt*coef/m);
}

int main(){
    Cuerpo Planeta;
    double t, dt=0.001;
    double r0=5, v0, m0=1;
    double omega, T;
    
    omega=sqrt(GM)*pow(r0,-1.5);
    v0= omega*r0;
    T=2*M_PI/omega;

    //Debug
    //cout<<T<<" "<<v0<<" "<<omega<<endl;
    //cout<<" "<<Zeta<<" "<<" "<<Chi<<" "<<" "<<Lambda<<" "<<Coeficiente1<<" "<<Coeficiente2<<endl;

    //------------( x0,y0,z0,Vx0,Vy0,Vz0, m0, R0);
    Planeta.Inicie( r0, 0, 0,  0, 0.75*v0,  0, m0, 0.5);



    for(t=0; t<T; t+=dt){
        cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<endl;
        // Mover por PEFRL
        Planeta.Mueva_r(dt,Zeta);

        Planeta.CalculeFuerza();
        Planeta.Mueva_V(dt, Coeficiente1);

        Planeta.Mueva_r(dt, Chi);

        Planeta.CalculeFuerza();
        Planeta.Mueva_V(dt, Lambda);

        Planeta.Mueva_r(dt, Coeficiente2);

        Planeta.CalculeFuerza();
        Planeta.Mueva_V(dt, Lambda);

        Planeta.Mueva_r(dt, Chi);

        Planeta.CalculeFuerza();
        Planeta.Mueva_V(dt, Coeficiente1);

        Planeta.Mueva_r(dt,Zeta);
    }
    return 0;
}