//Mi Primer Programa
#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;

int main(){
    vector3D a, b, c;

    a.load(1,2,3);
    cout<<a.norm()<<endl;

    return 0;
}