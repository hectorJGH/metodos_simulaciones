#include <iostream>
#include <cmath>
#include "vector.h"

using namespace std;

//Constantes globales
const double G=1;
const int N=2;

//Constantes Forest-Roth
const double Theta=1/(2-pow(2.0,1.0/3));
const double ThetaU2 = Theta/2;
const double UmThetaU2 = (1 - Theta)/2;
const double Um2Theta = 1-2*Theta;

class Cuerpo;
class Colisionador;

class Cuerpo{
private:
    vector3D r,rold, V, F;
    double m,R;
public:
    void Inicie(double x0,double y0,double z0,
                    double Vx0,double Vy0, double Vz0,
                    double m0,double R0);
    void SumeFuerza(vector3D F0);
    void BorreFuerza(void);
    void Arranque(double dt);
    void Mueva_r(double dt, double coef);
    void Mueva_V(double dt, double coef);
    double Getx(void){return r.x();}
    double Gety(void){return r.y();}
    double Getz(void){return r.z();}
    friend class Colisionador;
};

void Cuerpo::Inicie(double x0,double y0,double z0,
                    double Vx0,double Vy0, double Vz0,
                    double m0,double R0){
    r.load(x0,y0,z0); V.load(Vx0,Vy0,Vz0); m=m0; R=R0;
}

void Cuerpo::SumeFuerza(vector3D F0){
    F+=F0;
}
void Cuerpo::BorreFuerza(void){
    F.load(0,0,0);
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

//----------Clase Colisionador----------
class Colisionador{
private:
public:
    void CalculeFuerzas(Cuerpo * Planeta);
    void CalculeFuerzaEntre(Cuerpo & Planeta1, Cuerpo & Planeta);
    ///Se ponen referencias & para que el objeto cambie con la función
};

void Colisionador::CalculeFuerzas(Cuerpo * Planeta){
    int i, j;
    //Borrar fuerzas
    for(i=0;i<N; i++)
    Planeta[i].BorreFuerza();
    //Calcular las fuerzas entre todas las parejas
    for (i=0; i<N;i++)
        for(j=i+1;j<N;j++)
            CalculeFuerzaEntre(Planeta[i],Planeta[j]);
}

void Colisionador::CalculeFuerzaEntre(Cuerpo & Planeta1, Cuerpo & Planeta2){
    vector3D r21,n,F1; double d, F;
    r21=Planeta2.r-Planeta1.r; d = r21.norm();
    n= r21/d;
    F = G * Planeta1.m * Planeta2.m * pow(d, -1.5);
    F1=n*F; Planeta1.SumeFuerza(F1);Planeta1.SumeFuerza(F1*(-1));
}


int main(){
    Cuerpo Planeta[N];
    Colisionador Newton;
    double m0=10, m1=1, r=11;
    double M=m0+m1, x0 = -m1*r/M, x1= m0*r/M;
    double omega=sqrt(G*M/(r*r*r)), T= 2*M_PI/omega, V0= omega*x0, V1=omega*x1;
    double t, dt=0.1;
    int i;

    //----------( x0, y0, z0, Vx0, Vy0, Vz0, m0, R0);
    Planeta[0].Inicie( x0, 0, 0, 0, V0, 0 , m0, 1);
    Planeta[1].Inicie( x1, 0, 0, 0, V1, 0 , m1, 0.15);


    for(t=0; t<T; t+=dt){
        cout<<Planeta[1].Getx()<<" "<<Planeta[1].Gety()<<endl;
        // Mover por Forest-Roth
        for(i=0; i<N; i++) Planeta[i].Mueva_r(dt, ThetaU2);
        Newton.CalculeFuerzas(Planeta);

        for(i=0; i<N; i++) Planeta[i].Mueva_V(dt, Theta);

        for(i=0; i<N; i++) Planeta[i].Mueva_r(dt,UmThetaU2);
        Newton.CalculeFuerzas(Planeta);

        for(i=0; i<N; i++) Planeta[i].Mueva_V(dt, Um2Theta);

        for(i=0; i<N; i++) Planeta[i].Mueva_r(dt,UmThetaU2);
        Newton.CalculeFuerzas(Planeta);


        for(i=0; i<N; i++) Planeta[i].Mueva_V(dt, Theta);

        for(i=0; i<N; i++) Planeta[i].Mueva_r(dt, ThetaU2);
    }
    return 0;
}