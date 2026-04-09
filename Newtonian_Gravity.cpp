
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>

const int B = 4;
double T0 = 0;
double T1 = 5;

double m[B] = {2.2, 0.8, 0.9, 0.4};

double x_init[B] = {-0.50, -0.60, 0.50, 0.50};
double y_init[B] = {0.10, -0.20, 0.10, 0.40};
double vx_init[B] = {-0.84, 1.86, -0.44, 1.15};
double vy_init[B] = {0.65, 0.70, -1.50, -1.60};

void acceleration(double X[], double Y[], double AX[], double AY[]) {
    
    for (int i = 0; i < B; ++i) { AX[i] = 0.0; AY[i] = 0.0; }

    for (int i = 0; i < B; ++i) {
        for (int j = i + 1; j < B; ++j) {
            double dx = X[j] - X[i];
            double dy = Y[j] - Y[i];
            double r2 = dx*dx + dy*dy;
            double inv_r3 = 1.0 / (std::sqrt(r2) * r2);
            double f = inv_r3;

            AX[i] += m[j] * f * dx;  AY[i] += m[j] * f * dy;
            AX[j] -= m[i] * f * dx;  AY[j] -= m[i] * f * dy;
        }
    } 
}

void write_solution(int N, const char* path, double xf[], double yf[]) {
    double x[B]; 
    double y[B]; 
    double vx[B]; 
    double vy[B];
    double ax[B];
    double ay[B];

    for (int i = 0; i < B; ++i) {
        x[i] = x_init[i];
        y[i] = y_init[i];
        vx[i] = vx_init[i]; 
        vy[i] = vy_init[i];
    }

    double h = (T1 - T0) / double(N);
    double t = T0;

    std::ofstream out(path);
    out << std::setprecision(16);
    
    out << "t,x0,y0,x1,y1,x2,y2,x3,y3\n";

    
    out << t << "," << x[0] << "," << y[0] << "," << x[1] << "," << y[1] << "," 
    << x[2] << "," << y[2] << "," << x[3] << "," << y[3] << "\n";
    

    acceleration(x, y, ax, ay);

    for (int i = 0; i < N; ++i) {
        // v(n+1/2)
        for (int k = 0; k < B; ++k) {
            vx[k] += 0.5 * h * ax[k];
            vy[k] += 0.5 * h * ay[k];
        }

        // x(n+1)
        for (int k = 0; k < B; ++k) {
            x[k] += h * vx[k];
            y[k] += h * vy[k];
        }
        t += h;

        // a(n+1)
        acceleration(x, y, ax, ay);

        // v(n+1)
        for (int k = 0; k < B; ++k) {
            vx[k] += 0.5 * h * ax[k];
            vy[k] += 0.5 * h * ay[k];
        }

        
        out << t;
        for (int k = 0; k < B; ++k) out << "," << x[k] << "," << y[k];
        out << "\n";
    }

    out.close();

    for (int i = 0; i < B; ++i) {
        xf[i] = x[i];
        yf[i] = y[i];
    }
}

int main() {
    int N = 50000; 
    double xf[B], yf[B];

    write_solution(N, "traj.csv", xf, yf);

    std::cout << std::scientific << std::setprecision(8);
    
        std::cout << "Planet" << 1 << ": x = "<< xf[0] << " y = " << yf[0] << "\n"
                  << "Planet" << 2 << ": x = "<< xf[1] << " y = " << yf[1] << "\n"
                  << "Planet" << 3 << ": x = "<< xf[2] << " y = " << yf[2] << "\n"
                  << "Planet" << 4 << ": x = "<< xf[3] << " y = " << yf[3] << "\n";
    return 0;
}
