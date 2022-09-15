#include <iostream>
#include <cmath>

using namespace std;

const double g=9.8;

class Cuerpo;

class Cuerpo{
    private:
    double x,y,omegax,omegay,Fx,Fy,m,R;
    public:
    void Inicie(double x0,double y0,double omegax0,double omegay0,double m0,double R0);
    void CalculeTorque(void);
    void Muevase(double dt);
    double Getx(void){return x;}
    double Gety(void){return y;}
};

void Cuerpo::Inicie(double x0,double y0,double omegax0,double omegay0,double m0,double R0){
    x=x0; y=y0; omegax=omegax0; omegay=omegay0; m=m0; R=R0;
}

void Cuerpo::CalculeTorque(void){
    Fx=0; Fy=g*m;
}

void Cuerpo::Muevase(double dt){
    x+=omegax*dt; y+=omegay*dt;
    omegax += Fx/m * dt; omegay += Fy/m * dt;
}

int main(){
    Cuerpo Balon;
    double t, dt=0.1;
    //----------( x0, y0, omegax0, omegay0, m0, R0);
    Balon.Inicie( 0, 0, 40, 30, 0.453, 0.15);
    for(t=0; t<3; t+=dt){
        cout<<Balon.Getx()<<" "<<Balon.Gety()<<endl;
        Balon.CalculeTorque();
        Balon.Muevase(dt);
    }
    return 0;
}