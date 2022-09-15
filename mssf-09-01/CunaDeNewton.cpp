#include <iostream>
#include <cmath>
//#include "vector.h"

using namespace std;

//Constantes globales
const double G=980;
const int N=5;
const double K=1e7;

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
    double theta, omega, tau; 
    double m,R, l, I, x0;
public:
    void Inicie(double theta0, double omega0, double tau0, double m0,double R0, double l0, double x00);
    void SumeTorque(double tau);
    void BorreTorque(void);
    void Mueva_theta(double dt, double coef);
    void Mueva_omega(double dt, double coef);
    void Dibujese(void);
    double Getx(void){return x0 + l*sin(theta);}
    double Gety(void){return -l*cos(theta);}
    double Gettheta(void){return theta;}
    double GetI(void){return I;}
    friend class Colisionador;
};

void Cuerpo::Inicie(double theta0, double omega0, double tau0, 
    double m0,double R0, double l0, double x00){
    theta=theta0; omega =omega0; tau=tau0; 
    m=m0; R=R0; l=l0; I=m0*l0*l0; x0=x00;
}

void Cuerpo::SumeTorque(double tau0){
    tau+=tau0;
}
void Cuerpo::BorreTorque(void){
    tau=0;
}


void Cuerpo::Mueva_theta(double dt, double coef){
    theta+= omega*(dt*coef);
}
void Cuerpo::Mueva_omega(double dt, double coef){
    omega+= tau*(dt*coef/m);
}

void Cuerpo::Dibujese(void){
  cout<<" , "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety()<<"+"<<R<<"*sin(t) ,";
  cout<<x0<<"+"<<l/7<<"*t*sin("<<theta<<") , -"<<l/7<<"*t*cos("<<theta<<")";
}

//----------Clase Colisionador----------
class Colisionador{
private:
public:
    void CalculeTorques(Cuerpo * Pendulo);
    void CalculeTorqueEntre(Cuerpo & Pendulo1, Cuerpo & Pendulo2);
};
///Se ponen referencias & para que el objeto cambie con la función
// * Es para aplicar a los vectores


void Colisionador::CalculeTorques(Cuerpo * Pendulo){
    int i, j;
    //Borrar Torques
    for(i=0;i<N; i++){
    Pendulo[i].BorreTorque();
    double tau_g=-Pendulo[i].l*Pendulo[i].m*G*sin(Pendulo[i].theta);
    Pendulo[i].SumeTorque(tau_g);
    }
    //Calcular las Torques entre todas las parejas
    for (i=0; i<N;i++)
        for(j=i+1;j<N;j++)
            CalculeTorqueEntre(Pendulo[i],Pendulo[j]);
}

void Colisionador::CalculeTorqueEntre(Cuerpo & Pendulo1, Cuerpo & Pendulo2){
    double s=(Pendulo2.Getx()+Pendulo2.R)-(Pendulo1.Getx()-Pendulo1.R); double F=0;
    if(s>0) F=K*pow(s,1.5);
    Pendulo1.SumeTorque(F*Pendulo1.l); Pendulo2.SumeTorque(-F*Pendulo2.l);


    //vector3D r21,n,ds_dt; double d, tau;
    //r21=Pendulo2.r-Pendulo1.r; d = r21.norm();
    //n= r21/d;
    //tau = G * Pendulo1.m * Pendulo2.m * pow(d, -2.0);
    //ds_dt=n*tau; Pendulo1.SumeTorque(ds_dt);Pendulo2.SumeTorque(ds_dt*(-1));
}

//----------Funciones globales
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'Pendulos.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-12:12]"<<endl;
  cout<<"set yrange[-12:0]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 60"<<endl;
}

void InicieCuadro(void){
    cout<<"plot 0,0 ";
}

void TermineCuadro(void){
    cout<<endl;
}


int main(){
    //------Pendulo
    Cuerpo Pendulo[N];
    Colisionador Newton;
    double x0=4, L=10, R=2, m=10;
    double theta0=M_PI*0.2, omega0=0, m0=m, R0=R, l0=L;
    double theta1=0, omega1=0, m1=m, R1=R, l1=L;

    //------Dibujo
    int i;
    double T=sqrt(L/G);
    double t,tmax=5*T,dt=0.00001;
    double tdibujo,tcuadro=T/100;

    double xini=x0*(N-1)/2.0;
    //----------(theta0, omega0, tau, m0, R0, l0, x00);
    Pendulo[0].Inicie(theta0, omega0, 0, m0, R0, l0, xini);
    for(i=1; i<N; i++) Pendulo[i].Inicie(theta1, omega1, 0, m1, R1, l1, xini-i*x0);

    InicieAnimacion();

    for(t=0, tdibujo=0; t<tmax; t+=dt,tdibujo+=dt){
        //cout<<Pendulo[1].Getx()<<" "<<Pendulo[1].Gety()<<endl;

        //Dibujar
        if(tdibujo>tcuadro){
        InicieCuadro();
        for(i=0;i<N;i++) Pendulo[i].Dibujese();
        TermineCuadro();
        tdibujo=0;
        }

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