#include <iostream>
#include <cmath>

using namespace std;

double f(double t, double x){
    return x;
}

void UnPasoDeRungeKutta4(double & t, double & x, double dt){
    double ds, di, dx3, dx4;
    ds= dt* f(t,x);
    di = dt* f(t+dt/2,x+ds/2);
    dx3 = dt* f(t+dt/2,x+di/2);
    dx4 = dt* f(t+dt/2,x+dx3);

    x+= (ds + 2*di + 2*dx3 + dx4)/6;
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