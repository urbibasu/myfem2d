#include <cmath>
#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "rheology.hpp"
#include "utils.hpp"


static void elastic(double bulkm, double shearm, const double* de, double* s)
{
    /* increment the stress s according to the incremental strain de */
    double lambda = bulkm - 2. /3 * shearm;
    double dev = trace(de);

    for (int i=0; i<NDIMS; ++i)
        s[i] += 2 * shearm * de[i] + lambda * dev;
    for (int i=NDIMS; i<NSTR; ++i)
        s[i] += 2 * shearm * de[i];
}


void update_stress(const Variables& var, tensor_t& stress,
                   double_vec& stressyy,
                   tensor_t& strain, double_vec& plstrain,
                   double_vec& delta_plstrain, tensor_t& strain_rate)
{
    const int rheol_type = var.mat->rheol_type;

    for (int e=0; e<var.nelem; ++e) {
        // stress, strain and strain_rate of this element
        double* s = stress[e];
        double& syy = stressyy[e];
        double* es = strain[e];
        double* edot = strain_rate[e];

        // anti-mesh locking correction on strain rate
        if(1){
            double div = trace(edot);
            //double div2 = ((*var.volume)[e] / (*var.volume_old)[e] - 1) / var.dt;
            for (int i=0; i<NDIMS; ++i) {
                edot[i] += ((*var.edvoldt)[e] - div) / NDIMS;  // XXX: should NDIMS -> 3 in plane strain?
            }
        }

        // update strain with strain rate
        for (int i=0; i<NSTR; ++i) {
            es[i] += edot[i] * var.dt;
        }

        // modified strain increment
        double de[NSTR];
        for (int i=0; i<NSTR; ++i) {
            de[i] = edot[i] * var.dt;
        }

        switch (rheol_type) {
        case MatProps::rh_elastic:
            {
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);
                elastic(bulkm, shearm, de, s);
            }
            break;
        default:
            std::cerr << "Error: unknown rheology type: " << rheol_type << "\n";
            std::exit(1);
            break;
        }
        // std::cerr << "stress " << e << ": ";
        // print(std::cerr, s, NSTR);
        // std::cerr << '\n';
    }
}
