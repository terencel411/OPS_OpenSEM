void print_eddy_pos(double* x, double* y, double* z){
    for(int i{0}; i < eddies[0]; i++){
        printf("%f, %f, %f \n", x[i], y[i], z[i]);
    }
}

void print_eddy_info(double* x, double* y, double* z, double* r){
    for(int i{0}; i < eddies[0]; i++){
        printf("%f, %f, %f, %f \n", x[i], y[i], z[i], r[i]);
    }
}

void print_norm_x(double* x, double* r){
    for(int i{0}; i < eddies[0]; i++){
        printf("%f \n", x[i] / r[i]);
    }
}

void calc_eddies(int* eddies, const double* vol, const double* radii, const int& nregions){
    for(int i{0}; i < nregions; i++){
        //eddies[i] = trunc(pow(vol[i]/((x_max-x_min)*radii[3*i+1]*radii[3*i+2]), 1.0));
        eddies[i] = trunc(vol[i]/(radii[3*i]*radii[3*i+1]*radii[3*i+2]));
        //eddies[i] = trunc(pow(vol[i]/(radii[3*i]*radii[3*i+1]*radii[3*i+2]), 1.0));
        //eddies[i] = trunc(vol[i]/(x_radii[i]*radii[3*i+1]*radii[3*i+2]));
    }
}