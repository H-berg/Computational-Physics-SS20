#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

double a_direkt(double x){ //direkt berechnet
    return 1/sqrt(x)-1/sqrt(x+1);
}

double a_verb(double x){ //verbessert
    return 1/(sqrt(x+1)*sqrt(x))*1/(sqrt(x+1)+sqrt(x));
}

double b_direkt(double x){
    return (1-cos(x))*1/sin(x);
}

double b_verb(double x){
    return tan(x/2);
}

double c_direkt(double x, double y){
    return sin(x+y)-sin(x);
}

double c_verb(double x, double y){
    return 2*cos(x+y/2)*sin(y/2);
}

void a(){ //Aufgabenteil a
    std::ofstream file;
    file.open( "a_2_a.txt", std::ios_base::trunc );
    file << "#direkt \t verbessert \t relativer Fehler \n" ;
    double a_1, a_2;
    double rel_err_a; //relativer Fehler
    for(int i=pow(10,5); i<=pow(10, 6); i++){
        a_1 = a_direkt(i);
        a_2 = a_verb(i);
        rel_err_a = abs(a_1-a_2)/abs(a_2);
        file << a_1 << '\t' << a_2 << '\t' << rel_err_a << '\n';
    }
    file.flush();
    file.close();
}

void b(){ //Aufgabenteil b
    std::ofstream file;
    file.open( "a_2_b.txt", std::ios_base::trunc );
    file << "#direkt \t verbessert \t relativer Fehler \n" ;
    double b_1, b_2;
    double rel_err_b; //relativer Fehler
    double i = 0.000000001;
    while (i < 0.0000001){
        i = i+0.000000001;
        b_1 = b_direkt(i);
        b_2 = b_verb(i);
        rel_err_b = abs(b_1-b_2)/abs(b_2);
        file << b_1 << '\t' << b_2 << '\t' << rel_err_b << '\n';
    }
    file.flush();
    file.close();
}

void c(){ //Aufgabenteil c
    std::ofstream file;
    file.open( "a_2_c.txt", std::ios_base::trunc );
    file << "#direkt \t verbessert \t relativer Fehler \n" ;
    double c_1, c_2;
    double rel_err_c; //relativer Fehler
    double i = 1000000;
    double j = 0.000000000001; //Wert für x
    while (j <= 0.000000001){
        j = j+0.000000000001;
        c_1 = c_direkt(i, j);
        c_2 = c_verb(i, j);
        rel_err_c = abs(c_1-c_2)/abs(c_2);
        file << c_1 << '\t' << c_2 << '\t' << rel_err_c << '\n';
    }
    file.flush();
    file.close();
}

int main() {
    a();
    b();
    c();
    return 0;
}
