#include <iostream>
#include <cmath>

using namespace std;

double f(double t, double x){
    return x;
}

void UnPasoDeRungeKutta4(double & t, double & x, double dt){
    double dx1, dx2, dx3, dx4;
    dx1= dt* f(t,x);
    dx2 = dt* f(t+dt/2,x+dx1/2);
    dx3 = dt* f(t+dt/2,x+dx2/2);
    dx4 = dt* f(t+dt/2,x+dx3);

    x+= (dx1 + 2*dx2 + 2*dx3 + dx4)/6;
    t+=dt;
}

int main(){
    double t, x; double dt=0.01;

    for(t=0,x=1; t<2+dt/2; ){
        cout<<t<<" "<<x<<endl;
        UnPasoDeRungeKutta4(t,x,dt);
    }
    return 0;
}