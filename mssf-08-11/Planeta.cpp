#include <iostream>
#include <cmath>

using namespace std;

const double GM=1;

class Cuerpo;

class Cuerpo{
    private:
    double x,y,omegax,omegay,Fx,Fy,m,R,r;
    public:
    void Inicie(double x0,double y0,double omegax0,double omegay0,double m0,double R0, double r0);
    void CalculeTorque(void);
    void Muevase(double dt);
    double Getx(void){return x;}
    double Gety(void){return y;}
};

void Cuerpo::Inicie(double x0,double y0,double omegax0,double omegay0,double m0,double R0, double r0){
    x=x0; y=y0; omegax=omegax0; omegay=omegay0; m=m0; R=R0; r=r0;
}

void Cuerpo::CalculeTorque(void){
    double tau = -GM * m/(r*r);
    Fx=tau*x/r; Fy=tau*y/r;
}

void Cuerpo::Muevase(double dt){
    x+=omegax*dt; y+=omegay*dt;
    omegax += Fx/m * dt; omegay += Fy/m * dt;
}

int main(){
    Cuerpo Pendulo;
    double t, dt=0.1;
    double omega, r0=5, v0, T,m0=1;
    
    omega=sqrt(GM/(r0*r0*r0)); //en clase
    v0= omega*r0;
    T=2*M_PI/omega;

    //----------( x0, y0, omegax0, omegay0, m0, R0, r);
    Pendulo.Inicie( r0, 0,  0, v0, m0, 0.15, r0);
    
    for(t=0; t<T; t+=dt){
        cout<<Pendulo.Getx()<<" "<<Pendulo.Gety()<<endl;
        Pendulo.CalculeTorque();
        Pendulo.Muevase(dt);
    }
    return 0;
}