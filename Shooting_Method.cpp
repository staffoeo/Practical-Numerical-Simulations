
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

double a(double t, double x, double v) { return - (t * x * (t + 2.0)) / (2.0 + t*t*x*x); }

// RK4
double end_x(double s) {
    double t0 = 0.0; 
    double t1 = 10.0; 
    double x = 0.75; 
    double v = s;
    int N = 5000;
    double h = (t1 - t0) / N;
    double t = t0;
     
    for (int i = 0; i < N; ++i) {
        double k1x = h*v;
        double k1v = h*a(t, x, v);

        double x2 = x + 0.5*k1x;
        double v2 = v + 0.5*k1v;
        double k2x = h*v2;
        double k2v = h*a(t + 0.5*h, x2, v2);

        double x3 = x + 0.5*k2x;
        double v3 = v + 0.5*k2v;
        double k3x = h*v3;
        double k3v = h*a(t + 0.5*h, x3, v3);

        double x4 = x + k3x;
        double v4 = v + k3v;
        double k4x = h*v4;
        double k4v = h*a(t + h, x4, v4);

        x += (1/6.0)*(k1x + 2.0*k2x + 2.0*k3x + k4x);
        v += (1/6.0)*(k1v + 2.0*k2v + 2.0*k3v + k4v);
        t += h;
    }
    return x;
}

// trajectory for plots
void write_solution(double s, const char* path) {
    double t0 = 0.0;
    double t1 = 10.0; 
    double x = 0.75;
    double v = s;
    int N = 5000; 
    double h = (t1 - t0) / N;
    double t = t0;
    std::ofstream out(path);
    out << std::setprecision(16) << "t,x,v\n";
    out << t << "," << x << "," << v << "\n";
    for (int i = 0; i < N; ++i) {
        double k1x = h*v;
        double k1v = h*a(t, x, v);

        double x2 = x + 0.5*k1x;
        double v2 = v + 0.5*k1v;
        double k2x = h*v2;
        double k2v = h*a(t + 0.5*h, x2, v2);

        double x3 = x + 0.5*k2x;
        double v3 = v + 0.5*k2v;
        double k3x = h*v3;
        double k3v = h*a(t + 0.5*h, x3, v3);

        double x4 = x + k3x;
        double v4 = v + k3v;
        double k4x = h*v4;
        double k4v = h*a(t + h, x4, v4);

        x += (1/6.0)*(k1x + 2.0*k2x + 2.0*k3x + k4x);
        v += (1/6.0)*(k1v + 2.0*k2v + 2.0*k3v + k4v);
        t += h;
        out << t << "," << x << "," << v << "\n";
    }
}

// we want F(s) = 0
double F(double s) { return end_x(s) + 1.0; }


int main() {
    std::cout << std::scientific << std::setprecision(6);
    int k = 0;
    double step = 0.1;
    for (double sL = -4.0; sL<4.0; sL += step) {
        double sR = sL + step;
        double fL = F(sL);
        double fR = F(sR);
        if (fL * fR < 0.0) {
            for (int j = 0; j < 100; ++j) {
                double sM = 0.5*(sL + sR);
                double fM = F(sM);
                if (fL * fM < 0.0) {sR = sM; fR = fM;}
                else               {sL = sM; fL = fM;}

            }
        double s = 0.5 * (sR + sL);
        ++ k;
        std::cout << "s" << k << " = " << s << "\n";
        std::string fname = "solution" + std::to_string(k) + ".csv";
        write_solution(s, fname.c_str());
        }

    }
    return 0;
}

