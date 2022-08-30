#include <iostream>
#include <cmath>
#include "vector.h"

using namespace std;

const double GM=1;

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
    void Muevase(double dt);
    double Getx(void){return r.x();}
    double Gety(void){return r.y();}
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

void Cuerpo::Muevase(double dt){
    V += F*(dt/m);
    r+= V*dt;
}

int main(){
    Cuerpo Planeta;
    double t, dt=0.1;
    double omega, r0=5, v0, T,m0=1;
    
    omega=sqrt(GM/(r0*r0*r0)); //en clase
    v0= omega*r0;
    T=2*M_PI/omega;

    //----------( x0, y0, Vx0, Vy0, m0, R0);
    Planeta.Inicie( r0, 0, 0, 0, v0, 0 , m0, 0.15);
    Planeta.CalculeFuerza();
    Planeta.Arranque(dt);


    for(t=0; t<T; t+=dt){
        cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<endl;
        Planeta.CalculeFuerza();
        Planeta.Muevase(dt);
    }
    return 0;
}