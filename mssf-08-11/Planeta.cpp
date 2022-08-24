#include <iostream>
#include <cmath>

using namespace std;

const double GM=1;

class Cuerpo;

class Cuerpo{
    private:
    double x,y,Vx,Vy,Fx,Fy,m,R,r;
    public:
    void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0, double r0);
    void CalculeFuerza(void);
    void Muevase(double dt);
    double Getx(void){return x;}
    double Gety(void){return y;}
};

void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0, double r0){
    x=x0; y=y0; Vx=Vx0; Vy=Vy0; m=m0; R=R0; r=r0;
}

void Cuerpo::CalculeFuerza(void){
    double F = -GM * m/(r*r);
    Fx=F*x/r; Fy=F*y/r;
}

void Cuerpo::Muevase(double dt){
    x+=Vx*dt; y+=Vy*dt;
    Vx += Fx/m * dt; Vy += Fy/m * dt;
}

int main(){
    Cuerpo Planeta;
    double t, dt=0.1;
    double omega, r0=5, v0, T,m0=1;
    
    omega=sqrt(GM/(r0*r0*r0)); //en clase
    v0= omega*r0;
    T=2*M_PI/omega;

    //----------( x0, y0, Vx0, Vy0, m0, R0, r);
    Planeta.Inicie( r0, 0,  0, v0, m0, 0.15, r0);
    
    for(t=0; t<T; t+=dt){
        cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<endl;
        Planeta.CalculeFuerza();
        Planeta.Muevase(dt);
    }
    return 0;
}