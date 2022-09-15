#include <iostream>
#include <cmath>
#include "vector.h"

using namespace std;

const double GM=1;

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
    void Muevase(double dt);
    double Getx(void){return r.x();}
    double Gety(void){return r.y();}
};

void Cuerpo::Inicie(double x0,double y0,double z0,
                    double omegax0,double omegay0, double omegaz0,
                    double m0,double R0){
    r.load(x0,y0,z0); omega.load(omegax0,omegay0,omegaz0);
}

void Cuerpo::CalculeTorque(void){
    double normF = -GM * m/r.norm2();
    tau=r*normF/r.norm();
}

void Cuerpo::Arranque(double dt){
    rold=r;
    r=r + omega*dt + tau*(dt*dt/(2*m));
}

void Cuerpo::Muevase(double dt){
    vector3D rnew;
    rnew=2*r-rold + tau*(dt*dt/m); omega = (rnew - rold)/(2*dt); 
    rold=r;r=rnew;
}

int main(){
    Cuerpo Pendulo;
    double t, dt=0.1;
    double omega, r0=5, v0, T,m0=1;
    
    omega=sqrt(GM/(r0*r0*r0)); //en clase
    v0= omega*r0;
    T=2*M_PI/omega;

    //----------( x0, y0, omegax0, omegay0, m0, R0);
    Pendulo.Inicie( r0, 0, 0, 0, v0, 0 , m0, 0.15);
    Pendulo.CalculeTorque();
    Pendulo.Arranque(dt);


    for(t=0; t<T; t+=dt){
        cout<<Pendulo.Getx()<<" "<<Pendulo.Gety()<<endl;
        Pendulo.CalculeTorque();
        Pendulo.Muevase(dt);
    }
    return 0;
}