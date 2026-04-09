#include <cmath>
#include <cstdio>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

int idx(int i, int j, int nx) { return i + j*nx; }
double xcoord(int i, double hx) { return i*hx; }
double ycoord(int j, double hy) { return j*hy; }

void set_outer_boundary(std::vector<double>& phi, std::vector<int>& phi_fixed,
                        int nx, int ny, double hx, double hy) 
{
    
    for (int i = 0; i < nx; ++i){
        phi[idx(i, 0, nx)] = 0.0;
        phi_fixed[idx(i, 0, nx)] = 0;   
    }

    
    for (int j = 0; j < ny; ++j){
        phi[idx(0, j, nx)] = 0.0;
        phi_fixed[idx(0, j, nx)] = 0;   
    }

    
    for (int i = 0; i < nx; ++i){
        double x = xcoord(i, hx);
        phi[idx(i, ny-1, nx)] = std::sqrt(x);
        phi_fixed[idx(1, ny-1, nx)] = 1;   
    }

    
    for (int j = 0; j < ny; ++j){
        double y = ycoord(j, hy);
        phi[idx(nx-1, j, nx)] = std::sqrt(y);
        phi_fixed[idx(nx-1, j, nx)] = 1;
    }
}

void set_inner_boundary(std::vector<double>& phi, std::vector<int>& phi_fixed,
                        int nx, int ny, double hx, double hy) 
{
    const double AxL = 0.2, AxR = 0.6;
    const double AyB = 0.7, AyT = 0.9;
    const double tol1 = 0.5*std::min(hx,hy);

    for (int i = 0; i<nx; ++i) {
        double x = xcoord(i, hx);
        if (x < AxL - tol1 || x > AxR + tol1) continue;

        for (int j = 0; j<ny; ++j) {
            double y = ycoord(j, hy);
            if (y < AyB - tol1 || y > AyT + tol1) continue;

            bool x_inA_at_left   = (x >= AxL + tol1);
            bool x_inA_at_right  = (x <= AxR - tol1);
            bool y_inA_at_bottom = (y >= AyB + tol1);
            bool y_inA_at_top    = (y <= AyT - tol1);

            if (x_inA_at_left && x_inA_at_right && y_inA_at_bottom && y_inA_at_top) {
                phi[idx(i,j,nx)] = 0.0;
                phi_fixed[idx(i,j,nx)] = 1;
            }
        }
    }   

    for (int i = 0; i<nx; ++i) {
        double x = xcoord(i, hx);
        if (x < AxL - tol1 || x > AxR + tol1) continue;

        for (int j = 0; j<ny; ++j) {
            double y = ycoord(j, hy);
            if (y < AyB - tol1 || y > AyT + tol1) continue;

            bool x_onA_at_left   = std::abs(x-AxL) <= tol1 && (y >= AyB - tol1 && y <= AyT + tol1);
            bool x_onA_at_right  = std::abs(x-AxR) <= tol1 && (y >= AyB - tol1 && y <= AyT + tol1);
            bool y_onA_at_bottom = std::abs(y-AyB) <= tol1 && (x >= AxL - tol1 && x <= AxR + tol1);
            bool y_onA_at_top    = std::abs(y-AyT) <= tol1 && (x >= AxL - tol1 && x <= AxR + tol1);

            if (x_onA_at_left || x_onA_at_right || y_onA_at_bottom || y_onA_at_top) {
                phi[idx(i,j,nx)] = 1.0;
                phi_fixed[idx(i,j,nx)] = 1;
            }
        }
    }

    const double BxL = 0.6, BxR = 0.8;
    const double ByB = 0.1, ByT = 0.5;

    for (int j = 0; j<ny; ++j){
        double y = ycoord(j,hy);
        if ( y < ByB-tol1 || y > ByT+tol1) continue;

        for (int i=0; i<nx; ++i){
            double x = xcoord(i,hx);
            if (x < BxL - tol1 || x > BxR + tol1) continue;

            bool on_horizontal = std::fabs(y - ByT) <= tol1 && (x >= BxL - tol1 && x <= BxR + tol1);
            bool on_vertical   = std::fabs(x - BxR) <= tol1 && (y >= ByB - tol1 && y <= ByT + tol1);

            if (on_horizontal || on_vertical) {
                phi[idx(i,j,nx)] = -1.0;
                phi_fixed[idx(i,j,nx)] = 1;
            }
        }
    }
}

void sor_method(std::vector<double>& phi,
                const std::vector<int>& phi_fixed,
                int nx, int ny,
                double w, int max_iterations,
                double tol,
                int& it_used, double& final_change)
{
    it_used = 0;
    final_change = 0.0;

    for (int it = 1; it <= max_iterations; ++it) {
        double max_change = 0.0;

        
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                int id = idx(i, j, nx);
                if (phi_fixed[id] != 0) continue; 

                double right  = phi[idx(i+1, j, nx)];
                double left   = phi[idx(i-1, j, nx)];
                double top    = phi[idx(i, j+1, nx)];
                double bottom = phi[idx(i, j-1, nx)];

                double gs = 0.25 * (right + left + top + bottom);

                double old_val = phi[id];
                double new_val = (1.0 - w) * old_val + w * gs;

                double change = std::fabs(new_val - old_val);
                if (change > max_change) max_change = change;

                phi[id] = new_val;
            }
        }

        
        for (int j = 0; j < ny; ++j) {
            phi[idx(0, j, nx)] = phi[idx(1, j, nx)];
        }

        
        for (int i = 0; i < nx; ++i) {
            phi[idx(i, 0, nx)] = phi[idx(i, 1, nx)];
        }

        it_used = it;
        final_change = max_change;

        if (max_change < tol) {
            break;
        }
    }
}

int main() {
    
    const int nx = 101;
    const int ny = 101;

    const double tol = 1e-6;    
    const int    max_iterations = 3000;

    const double hx = 1.0 / (nx - 1);
    const double hy = 1.0 / (ny - 1);

    const int wsteps = 100;

    double w = 1.99;
        
    std::vector<double> phi(nx * ny, 0.0); 
    std::vector<int>    phi_fixed(nx * ny, 0); 

    set_outer_boundary(phi, phi_fixed, nx, ny, hx, hy);
    set_inner_boundary(phi, phi_fixed, nx, ny, hx, hy);

    int it_used = 0;
    double final_change = 0.0;
        
    sor_method(phi, phi_fixed, nx, ny, w, max_iterations, tol, it_used, final_change);

    int i0 = (int)std::lround(0.3 / hx);
    int j0 = (int)std::lround(0.5 / hy);
    double dphidy = (phi[idx(i0, j0+1, nx)] - phi[idx(i0, j0-1, nx)]) / (2.0 * hy); 

    std::cout << std::scientific << std::setprecision(6) << " dphi/dy(0.3,0.5) = " << dphidy << std::endl;

    return 0;
    
       
}