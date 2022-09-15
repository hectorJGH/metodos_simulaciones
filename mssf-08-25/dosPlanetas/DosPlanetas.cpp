#include <iostream>
#include <cmath>
#include "vector.h"

using namespace std;

//Constantes globales
const double G=1;
const int N=2;

//Constantes PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e00;
const double Chi=-0.6626458266981849e-1;

const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

class Cuerpo;
class Colisionador;

class Cuerpo{
private:
    vector3D r,rold, omega, tau;
    double m,R;
public:
    void Inicie(double x0,double y0,double z0,
                    double omegax0,double omegay0, double omegaz0,
                    double m0,double R0);
    void SumeTorque(vector3D F0);
    void BorreTorque(void);
    void Arranque(double dt);
    void Mueva_theta(double dt, double coef);
    void Mueva_omega(double dt, double coef);
    double Getx(void){return r.x();}
    double Gety(void){return r.y();}
    double Getz(void){return r.z();}
    friend class Colisionador;
};

void Cuerpo::Inicie(double x0,double y0,double z0,
                    double omegax0,double omegay0, double omegaz0,
                    double m0,double R0){
    r.load(x0,y0,z0); omega.load(omegax0,omegay0,omegaz0); m=m0; R=R0;
}

void Cuerpo::SumeTorque(vector3D F0){
    tau+=F0;
}
void Cuerpo::BorreTorque(void){
    tau.load(0.0,0.0,0.0);
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

//----------Clase Colisionador----------
class Colisionador{
private:
public:
    void CalculeTorques(Cuerpo * Pendulo);
    void CalculeTorqueEntre(Cuerpo & Pendulo1, Cuerpo & Pendulo2);
    ///Se ponen referencias & para que el objeto cambie con la función
    // * Es para aplicar a los vectores
};

void Colisionador::CalculeTorques(Cuerpo * Pendulo){
    int i, j;
    //Borrar Torques
    for(i=0;i<N; i++)
    Pendulo[i].BorreTorque();
    //Calcular las Torques entre todas las parejas
    for (i=0; i<N;i++)
        for(j=i+1;j<N;j++)
            CalculeTorqueEntre(Pendulo[i],Pendulo[j]);
}

void Colisionador::CalculeTorqueEntre(Cuerpo & Pendulo1, Cuerpo & Pendulo2){
    vector3D r21,n,ds_dt; double d, tau;
    r21=Pendulo2.r-Pendulo1.r; d = r21.norm();
    n= r21/d;
    tau = G * Pendulo1.m * Pendulo2.m * pow(d, -2.0);
    ds_dt=n*tau; Pendulo1.SumeTorque(ds_dt);Pendulo2.SumeTorque(ds_dt*(-1));
}


int main(){
    Cuerpo Pendulo[N];
    Colisionador Newton;
    double m0=10, m1=1, r=11;
    double M=m0+m1, x0 = -m1*r/M, s= m0*r/M;
    double omega=sqrt(G*M/(r*r*r)), T= 2*M_PI/omega, omega0= omega*x0, omega1=omega*s;
    double t, dt=0.1;
    int i;

    //----------( x0, y0, z0, omegax0, omegay0, omegaz0, m0, R0);
    Pendulo[0].Inicie( x0, 0, 0, 0, omega0, 0 , m0, 1);
    Pendulo[1].Inicie( s, 0, 0, 0, omega1, 0 , m1, 0.15);


    for(t=0; t<T; t+=dt){
        cout<<Pendulo[1].Getx()<<" "<<Pendulo[1].Gety()<<endl;
        // Mover por Forest-Roth
        for(i=0; i<N; i++) Pendulo[i].Mueva_theta(dt,Zeta);

        Newton.CalculeTorques(Pendulo);
        for(i=0; i<N; i++) Pendulo[i].Mueva_omega(dt, Coeficiente1);

        for(i=0; i<N; i++) Pendulo[i].Mueva_theta(dt,Chi);

        Newton.CalculeTorques(Pendulo);
        for(i=0; i<N; i++) Pendulo[i].Mueva_omega(dt, Lambda);

        for(i=0; i<N; i++) Pendulo[i].Mueva_theta(dt,Coeficiente2);

        Newton.CalculeTorques(Pendulo);
        for(i=0; i<N; i++) Pendulo[i].Mueva_omega(dt, Lambda);

        for(i=0; i<N; i++) Pendulo[i].Mueva_theta(dt,Chi);

        Newton.CalculeTorques(Pendulo);
        for(i=0; i<N; i++) Pendulo[i].Mueva_omega(dt, Coeficiente1);

        for(i=0; i<N; i++) Pendulo[i].Mueva_theta(dt,Zeta);
    }
    return 0;
}