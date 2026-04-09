#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>

double f(double t, double x) {
    return (t - 1.0) * (t - 1.0) * (x + 1.0);
}

double x_exact(double t) {
    return std::exp(t*t*t/3.0 - t*t + t) - 1.0; // exact solution
}

double rk4_N_steps(int N, double T) {
    double h = T / N;
    double t = 0.0;
    double x = 0.0;               // x(0) = 0
    int steps = 0;

    while (t < T) {
        double dt = std::min(h, T - t);           // safe last partial step if needed
        double k1 = f(t, x);
        double k2 = f(t + 0.5*dt, x + 0.5*dt*k1);
        double k3 = f(t + 0.5*dt, x + 0.5*dt*k2);
        double k4 = f(t + dt,     x + dt*k3);
        x += (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4);
        t += dt;
        ++steps;
        if (steps >= N && T - t <= 1e-15) break;  
    }
    return x; // x(T)
}

int main() {
    const double T = 2.0;          
    const int N_min = 1;
    const int N_max = 500;

    std::ofstream out("rk4_N_table.csv");
    out << std::fixed << std::setprecision(12);
    out << "N,h,x_N,x_exact_T,abs_error\n";

    const double xT = x_exact(T);

    std::cout << "N,h,x_N,abs_error\n";
    std::cout << std::scientific << std::setprecision(6);

    for (int N = N_min; N <= N_max; N+=1) {
        double h   = T / N;
        double xN  = rk4_N_steps(N, T);
        double err = std::fabs(xN - xT);

        out << N << "," << std::setprecision(12) << h << ","
            << xN << "," << xT << "," << err << "\n";

        if (N <= 12 || N % 50 == 0 || N == N_max) {
            std::cout << N << "," << h << "," << xN << "," << err << "\n";
        }
    }
    out.close();

   
    const double h_traj = 0.2;                 
    std::ofstream out2("rk4_solution.csv");
    out2 << std::fixed << std::setprecision(12);
    out2 << "t,x_exact,x_rk4,abs_error\n";

    {
        double t = 0.0, x = 0.0;
        out2 << t << "," << x_exact(t) << "," << x << "," << std::fabs(x - x_exact(t)) << "\n";

        while (t < T) {
            double dt = std::min(h_traj, T - t);  // clamp final step to avoid overshoot
            double k1 = f(t, x);
            double k2 = f(t + 0.5*dt, x + 0.5*dt*k1);
            double k3 = f(t + 0.5*dt, x + 0.5*dt*k2);
            double k4 = f(t + dt,     x + dt*k3);
            x += (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4);
            t += dt;

            double xe = x_exact(t);
            out2 << t << "," << xe << "," << x << "," << std::fabs(x - xe) << "\n";
        }
    }
    out2.close();

    return 0;
}


