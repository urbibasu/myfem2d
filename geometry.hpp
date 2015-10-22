#ifndef DYNEARTHSOL3D_GEOMETRY_HPP
#define DYNEARTHSOL3D_GEOMETRY_HPP

double dist2(const double* a, const double* b);
void compute_volume(const array_t &coord, const conn_t &connectivity,
                    double_vec &volume);

double compute_dt(const Param& param, const Variables& var);

void compute_shape_fn(const array_t &coord, const conn_t &connectivity,
                      const double_vec &volume,
                      shapefn &shpdx, shapefn &shpdy);

#endif
