#include <vector>
#include <algorithm>
#include <execution>
#include <iostream>
#include <array>
#include <cmath>
#include <string>
#include <ranges>
#include <fstream>
#include <iomanip>
#include <cassert>

using namespace std;

constexpr auto policy = execution::par;


int main() {
  //Dimensiones del sistema
  constexpr auto sigma_2  = 1.;
  constexpr auto sigma_1  = 2. * sigma_2;
  constexpr auto sigma_12 = ( sigma_1+sigma_2 ) / 2.;

  
  auto r = vector<double>(100);  //Cree un vector de floats de 101 zeros 
//  auto r = vector<double>{5.5}; 
  iota( begin(r), end(r), sigma_1*100. ); //Convertir el vector en {1000,1001,...,1099}
					       

  transform(policy, begin(r), end(r), begin(r), [](auto i){
		  return i/100;}); //Transformar en el vector {10.00,10.01,...,11.00}
 
  ofstream fout("beta_u_dis_bin_bi_s1_2.csv"); // Create/output to a file

  //Funcion a integrar 
  auto f = [sigma_12](auto rho, auto q){
  	return ((cyl_bessel_j(0, rho*q ) * cyl_bessel_j(1, q*sigma_12 )  * cyl_bessel_j(1, q*sigma_12 ) ) / q) ;
  };
  
  // Propiedades de la integración
  constexpr auto a         = 0.0001 ;
  constexpr auto b         = 10000.0001;
  constexpr auto h         = 0.0001;
  constexpr auto mitad_n = (b-a)/(2*h);

  for(auto rho:r) {
    // Inicializacion de indices pares e impares para aplicar el método de Simpson Compuesto
    auto index_par   = vector<double>(mitad_n - 1);
    ranges::iota( index_par, 0 );
    transform(policy, begin(index_par), end(index_par), begin(index_par), [h,a](auto i){return (a + h*(i*2+2));  });

    auto index_inpar   = vector<double>(mitad_n);
    ranges::iota( index_inpar, 0 );
    transform(policy, begin(index_inpar), end(index_inpar), begin(index_inpar), [h,a](auto i){return (a + h*(i*2+1));  });


    // Primer termino para la aproximación
    auto XI0 = f(rho, a) + f(rho,b);
    auto XI1 = transform_reduce(policy, begin(index_inpar), end(index_inpar), 0., plus{}, [f,rho](auto X){
		    return f(rho, X);
		    });
    auto XI2 = transform_reduce(policy, begin(index_par), end(index_par), 0., plus{}, [f,rho](auto X){
		    return f(rho, X);
		    });
    auto XI  = h*(XI0 + 2*XI2 + 4*XI1)/3;  

    cout << XI << endl;
    fout << setprecision(10) << rho << "," <<XI << endl; 
  }

  fout.close();
  return 0;
}
