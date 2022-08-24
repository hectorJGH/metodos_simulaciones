#include <iostream>
#include <cmath>

using namespace std;

double f(double x){
  return sin(x)/x;
}

const double ErrMax=1e-7;

int main(){
  double x;
  double a=2, b=4, m, fa, fm;

  fa=f(a);
  while(b-a>ErrMax){
    m = (b+a)/2; fm=f(m);
      if(fa*fm>0)
	{a=m; fa=fm;} //correr a hasta m;
      else
	b=m; //correr b hasta m;
    }
  cout<<"El cero es "<<(a+b)/2<<endl;
  
  return 0;
}
