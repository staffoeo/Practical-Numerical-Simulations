#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <array>
#include <filesystem>

using State = std::array<double, 4>;

State f(double x, const State& u){

    State du;
    du[0] = u[1]; // u1'=u2
    du[1] = u[2]; // u2'=u3
    du[2] = u[3]; // u3'=u4
    du[3] = -5*u[2] - 0.1*x*std::pow(u[0],3); // u4'=-5u3-0.1*x*u1^3
    return du;
}

int main() { 
    const double x0 = 0.0, X = 30.0; 
    const double h = 0.8; 
    double x = x0; 
    State u{1.0, 0.0, 0.0, 0.0}; // initial conditions: u1(0)=1, u2(0)=0, u3(0)=0, u4(0)=0
    
    std::ofstream file("pns3_solution.csv"); 
    file << std::fixed << std::setprecision(12); 
    file << "x,u1\n"; 
    file << x << "," << u[0] << "\n" ;

    while (x < X) { // loop until we reach X 
        double dx = h; 
        if (x + dx > X) dx = X - x; // last partial step if needed 
        State k1 = f(x, u); 
        State k2 = f(x + dx/2.0, {u[0] + dx*k1[0]/2.0, u[1] + dx*k1[1]/2.0, u[2] + dx*k1[2]/2.0, u[3] + dx*k1[3]/2.0}); 
        State k3 = f(x + dx/2.0, {u[0] + dx*k2[0]/2.0, u[1] + dx*k2[1]/2.0, u[2] + dx*k2[2]/2.0, u[3] + dx*k2[3]/2.0}); 
        State k4 = f(x + dx, {u[0] + dx*k3[0], u[1] + dx*k3[1], u[2] + dx*k3[2], u[3] + dx*k3[3]}); 
        
        for (int i=0; i<4; ++i) {
            u[i] += (dx/6.0)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]); // update u using RK4 formula component-wise 
            } 
        
        x += dx; // increment x 
            
        file << x << "," << u[0] << "\n"; 
    } 
                        
    file.close(); 
    std::cout << std::defaultfloat << std::setprecision(8); // typical default
    std::cout << x << " " << u[0] << " " << u[1] << " " << u[2] << " " << u[3] << "\n"; // printing these values now as the loop has completed fully and we are now at x = X = 30 so all of these values printed out are evaluated at x = X =30

    return 0; 
}
