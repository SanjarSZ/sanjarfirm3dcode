#include "regular_grid_interpolant_3d.h"
#include <xtensor/xarray.hpp>
#include "xtensor/xlayout.hpp"
#define _USE_MATH_DEFINES
#include <math.h>
#include <boost/format.hpp>

#define _EPS_ 1e-13

template<class Array>
const int RegularGridInterpolant3D<Array>::simdcount;

template<class Array>
void RegularGridInterpolant3D<Array>::interpolate_batch(std::function<Vec(Vec, Vec, Vec)> &f) {
    std::cout << "DEBUG: interpolate_batch() called - this should NOT happen during load!" << std::endl;
    std::cout << "DEBUG: This is the expensive computation that should be avoided!" << std::endl;
    int BATCH_SIZE = 16384;
    int NUM_BATCHES = dofs_to_keep/BATCH_SIZE + (dofs_to_keep % BATCH_SIZE != 0);
    for (int i = 0; i < NUM_BATCHES; ++i) {
        uint32_t first = i * BATCH_SIZE;
        uint32_t last = std::min((uint32_t)((i+1) * BATCH_SIZE), dofs_to_keep);
        Vec xsub(xdoftensor_reduced.begin() + first, xdoftensor_reduced.begin() + last);
        Vec ysub(ydoftensor_reduced.begin() + first, ydoftensor_reduced.begin() + last);
        Vec zsub(zdoftensor_reduced.begin() + first, zdoftensor_reduced.begin() + last);
        Vec fxyzsub  = f(xsub, ysub, zsub);
        for (int j = 0; j < last-first; ++j) {
            for (int l = 0; l < value_size; ++l) {
                vals[first * value_size + j * value_size + l] = fxyzsub[j * value_size + l];
            }
        }
    }
    int degree = rule.degree;
    all_local_vals_map = std::unordered_map<int, AlignedPaddedVec>();
    all_local_vals_map.reserve(cells_to_keep);

    for (int xidx = 0; xidx < nx; ++xidx) {
        for (int yidx = 0; yidx < ny; ++yidx) {
            for (int zidx = 0; zidx < nz; ++zidx) {
                int meshidx = idx_cell(xidx, yidx, zidx);
                if(skip_cell[meshidx])
                    continue;
                AlignedPaddedVec local_vals(local_vals_size, 0.);
                for (int i = 0; i < degree+1; ++i) {
                    for (int j = 0; j < degree+1; ++j) {
                        for (int k = 0; k < degree+1; ++k) {
                            int offset = value_size*full_to_reduced_map[idx_dof(xidx*degree+i, yidx*degree+j, zidx*degree+k)];
                            int offset_local = padded_value_size * idx_dof_local(i, j, k);
                            for (int l = 0; l < value_size; ++l) {
                                local_vals[offset_local + l] = vals[offset + l];
                            }
                        }
                    }
                }
                all_local_vals_map.insert({meshidx, local_vals});
            }
        }
    }
}

template<class Array>
void RegularGridInterpolant3D<Array>::evaluate_batch(Array& xyz, Array& fxyz){
    if(fxyz.layout() != xt::layout_type::row_major)
          throw std::runtime_error("fxyz needs to be in row-major storage order");
    int npoints = xyz.shape(0);
    for (int i = 0; i < npoints; ++i) {
        evaluate_inplace(xyz(i, 0), xyz(i, 1), xyz(i, 2), fxyz.data() + value_size*i);
    }
}

template<class Array>
void RegularGridInterpolant3D<Array>::evaluate_batch_1D(Array& xyz, Array& fxyz){
    if(fxyz.layout() != xt::layout_type::row_major)
          throw std::runtime_error("fxyz needs to be in row-major storage order");
    int npoints = xyz.shape(0);
    for (int i = 0; i < npoints; ++i) {
        evaluate_inplace(xyz(i, 0), fxyz.data() + value_size*i);
    }
}

template<class Array>
Vec RegularGridInterpolant3D<Array>::evaluate(double x, double y, double z){
    Vec fxyz(value_size, 0.);
    evaluate_inplace(x, y, z, fxyz.data());
    return fxyz;

}

template<class Array>
int RegularGridInterpolant3D<Array>::locate_unsafe(double x, double y, double z){
    int xidx = int(nx*(x-xmin)/(xmax-xmin)); // find idx so that xmesh[xidx] <= x <= xs[xidx+1]
    int yidx = int(ny*(y-ymin)/(ymax-ymin));
    int zidx = int(nz*(z-zmin)/(zmax-zmin));
    return idx_cell(xidx, yidx, zidx);
}

template<class Array>
void RegularGridInterpolant3D<Array>::evaluate_inplace(double x, double y, double z, double* res){

    // to avoid funny business when the data is just a tiny bit out of bounds
    // due to machine precision, we perform this check and shift
    if(x >= xmax) x -= _EPS_;
    else if (x <= xmin) x += _EPS_;
    if(y >= ymax) y -= _EPS_;
    else if (y <= ymin) y += _EPS_;
    if(z >= zmax) z -= _EPS_;
    else if (z <= zmin) z += _EPS_;

    int xidx = int(nx*(x-xmin)/(xmax-xmin)); // find idx so that xmesh[xidx] <= x <= xs[xidx+1]
    int yidx = int(ny*(y-ymin)/(ymax-ymin));
    int zidx = int(nz*(z-zmin)/(zmax-zmin));
    if(!out_of_bounds_ok){
        if(xidx < 0 || xidx >= nx)
            throw std::runtime_error((boost::format("xidxs={} not within [0, {}]") % xidx % (nx-1)).str());
        if(yidx < 0 || yidx >= ny)
            throw std::runtime_error((boost::format("yidxs={} not within [0, {}]") % yidx % (ny-1)).str());
        if(zidx < 0 || zidx >= nz)
            throw std::runtime_error((boost::format("zidxs={} not within [0, {}]") % zidx % (nz-1)).str());
    }
    double xlocal = (x-xmesh[xidx])/hx;
    double ylocal = (y-ymesh[yidx])/hy;
    double zlocal = (z-zmesh[zidx])/hz;
    return evaluate_local(xlocal, ylocal, zlocal, idx_cell(xidx, yidx, zidx), res);
}

template<class Array>
void RegularGridInterpolant3D<Array>::evaluate_inplace(double x, double* res){

    // to avoid funny business when the data is just a tiny bit out of bounds
    // due to machine precision, we perform this check and shift
    if(x >= xmax) x -= _EPS_;
    else if (x <= xmin) x += _EPS_;

    int xidx = int(nx*(x-xmin)/(xmax-xmin)); // find idx so that xmesh[xidx] <= x <= xs[xidx+1]

    if(!out_of_bounds_ok){
        if(xidx < 0 || xidx >= nx)
            throw std::runtime_error((boost::format("xidxs={} not within [0, {}]") % xidx % (nx-1)).str());
    }
    double xlocal = (x-xmesh[xidx])/hx;
    return evaluate_local(xlocal, idx_cell(xidx, 0, 0), res);
}

template<class Array>
void RegularGridInterpolant3D<Array>::evaluate_local(double x, double y, double z, int cell_idx, double* res)
{
    int degree = rule.degree;
    auto got = all_local_vals_map.find(cell_idx);
    if (got == all_local_vals_map.end()) {
        if(out_of_bounds_ok)
            return;
        else
            throw std::runtime_error((boost::format("cell_idx={} not in all_local_vals_map") % cell_idx).str());
    }

    double* vals_local = got->second.data();
    if(xsimd::simd_type<double>::size >= 3){
        simd_t xyz;
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;
        for (int k = 0; k < degree+1; ++k) {
            simd_t temp = this->rule.basis_fun(k, xyz);
            pkxs[k] = temp[0];
            pkys[k] = temp[1];
            pkzs[k] = temp[2];
        }
    } else {
        for (int k = 0; k < degree+1; ++k) {
            pkxs[k] = this->rule.basis_fun(k, x);
            pkys[k] = this->rule.basis_fun(k, y);
            pkzs[k] = this->rule.basis_fun(k, z);
        }
    }

    // Potential optimization: use barycentric interpolation here right now the
    // implementation in O(degree^3) in memory and O(degree^4) in computation,
    // using Barycentric interpolation this could be reduced to O(degree^3) in
    // memory and O(degree^3) in computation.
    for(int l=0; l<padded_value_size; l += simdcount) {
        simd_t sumi(0.);
        int offset_local = l;
        double* val_ptr = &(vals_local[offset_local]);
        for (int i = 0; i < degree+1; ++i) {
            simd_t sumj(0.);
            for (int j = 0; j < degree+1; ++j) {
                simd_t sumk(0.);
                for (int k = 0; k < degree+1; ++k) {
                    double pkz = pkzs[k];
                    sumk = xsimd::fma(xsimd::load_aligned(val_ptr), simd_t(pkz), sumk);
                    val_ptr += padded_value_size;
                }
                double pjy = pkys[j];
                sumj = xsimd::fma(sumk, simd_t(pjy), sumj);
            }
            double pix = pkxs[i];
            sumi = xsimd::fma(sumj, simd_t(pix), sumi);
        }
        for (int ll = 0; ll < std::min(simdcount, value_size-l); ++ll) {
            res[l+ll] = sumi[ll];
        }
    }
}
//TODO memory usage not fixed
template<class Array>
void RegularGridInterpolant3D<Array>::evaluate_local(double x, int cell_idx, double* res)
{
    int degree = rule.degree;
    auto got = all_local_vals_map.find(cell_idx);
    if (got == all_local_vals_map.end()) {
        if(out_of_bounds_ok)
            return;
        else
            throw std::runtime_error((boost::format("cell_idx={} not in all_local_vals_map") % cell_idx).str());
    }

    double* vals_local = got->second.data();

    for (int k = 0; k < degree+1; ++k) {
        pkxs[k] = this->rule.basis_fun(k, x);
    }

    for(int l=0; l<padded_value_size; l += simdcount) {
        simd_t sumi(0.);
        int offset_local = l;
        double* val_ptr = &(vals_local[offset_local]);
        for (int i = 0; i < degree+1; ++i) {
            double pkx = pkxs[i];
            sumi = xsimd::fma(xsimd::load_aligned(val_ptr), simd_t(pkx), sumi);
            val_ptr += padded_value_size * (degree+1) * (degree+1);
        }
        for (int ll = 0; ll < std::min(simdcount, value_size-l); ++ll) {
            res[l+ll] = sumi[ll];
        }
    }
}

template<class Array>
std::pair<double, double> RegularGridInterpolant3D<Array>::estimate_error(std::function<Vec(Vec, Vec, Vec)> &f, int samples) {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, +1.0);
    double err = 0;
    double errsq = 0;
    Vec xs(samples, 0.);
    Vec ys(samples, 0.);
    Vec zs(samples, 0.);
    Array xyz = xt::zeros<double>({samples, 3});
    Array fhxyz = xt::zeros<double>({samples, value_size});
    for (int i = 0; i < samples; ++i) {
        xs[i] = xmin + distribution(generator)*(xmax-xmin);
        ys[i] = ymin + distribution(generator)*(ymax-ymin);
        zs[i] = zmin + distribution(generator)*(zmax-zmin);
        xyz(i, 0) = xs[i];
        xyz(i, 1) = ys[i];
        xyz(i, 2) = zs[i];
    }
    Vec fx = f(xs, ys, zs);
    this->evaluate_batch(xyz, fhxyz);
    for (int i = 0; i < samples; ++i) {
        double diff = 0.;
        for (int l = 0; l < value_size; ++l) {
            diff += std::pow(fx[value_size*i+l]-fhxyz(i, l), 2);
        }
        diff = std::sqrt(diff);
        err += diff;
        errsq += diff*diff;
    }
    double mean = err/samples;
    double std = std::sqrt((errsq - err*err/samples)/(samples-1)/samples);
    return std::make_pair(mean-std, mean+std);
}



Vec linspace(double min, double max, int n, bool endpoint) {
    Vec res(n, 0.);
    if(endpoint) {
        double h = (max-min)/(n-1);
        for (int i = 0; i < n; ++i)
            res[i] = min + i*h;
    } else {
        double h = (max-min)/n;
        for (int i = 0; i < n; ++i)
            res[i] = min + i*h;
    }
    return res;
}


template<class Array>
std::map<std::string, std::vector<double>> RegularGridInterpolant3D<Array>::get_interpolant_data() const {
    std::map<std::string, std::vector<double>> data;
    
    // Save the interpolated values
    data["vals"] = vals;
    
    // Save grid information
    data["xmesh"] = xmesh;
    data["ymesh"] = ymesh;
    data["zmesh"] = zmesh;
    data["xdof"] = xdof;
    data["ydof"] = ydof;
    data["zdof"] = zdof;
    data["xdoftensor_reduced"] = xdoftensor_reduced;
    data["ydoftensor_reduced"] = ydoftensor_reduced;
    data["zdoftensor_reduced"] = zdoftensor_reduced;
    
    // Save grid parameters
    data["nx"] = {static_cast<double>(nx)};
    data["ny"] = {static_cast<double>(ny)};
    data["nz"] = {static_cast<double>(nz)};
    data["hx"] = {hx};
    data["hy"] = {hy};
    data["hz"] = {hz};
    data["xmin"] = {xmin};
    data["ymin"] = {ymin};
    data["zmin"] = {zmin};
    data["xmax"] = {xmax};
    data["ymax"] = {ymax};
    data["zmax"] = {zmax};
    data["value_size"] = {static_cast<double>(value_size)};
    data["padded_value_size"] = {static_cast<double>(padded_value_size)};
    data["dofs_to_keep"] = {static_cast<double>(dofs_to_keep)};
    data["cells_to_keep"] = {static_cast<double>(cells_to_keep)};
    data["local_vals_size"] = {static_cast<double>(local_vals_size)};
    
    // Save interpolation rule
    data["rule_degree"] = {static_cast<double>(rule.degree)};
    data["rule_nodes"] = rule.nodes;
    data["rule_scalings"] = rule.scalings;
    
    // Save mappings
    data["reduced_to_full_map"] = std::vector<double>(reduced_to_full_map.begin(), reduced_to_full_map.end());
    data["full_to_reduced_map"] = std::vector<double>(full_to_reduced_map.begin(), full_to_reduced_map.end());
    
    // Save skip information
    data["skip_cell"] = std::vector<double>(skip_cell.begin(), skip_cell.end());
    
    // Save local values if available
    if (!all_local_vals_map.empty()) {
        for (const auto& pair : all_local_vals_map) {
            std::string key = "local_vals_" + std::to_string(pair.first);
            data[key] = std::vector<double>(pair.second.begin(), pair.second.end());
        }
    }
    
    return data;
}

template<class Array>
void RegularGridInterpolant3D<Array>::set_interpolant_data(const std::map<std::string, std::vector<double>>& data) {
    // Load the interpolated values
    if (data.find("vals") != data.end()) {
        vals = data.at("vals");
    }
    
    // Load grid information
    if (data.find("xmesh") != data.end()) xmesh = data.at("xmesh");
    if (data.find("ymesh") != data.end()) ymesh = data.at("ymesh");
    if (data.find("zmesh") != data.end()) zmesh = data.at("zmesh");
    if (data.find("xdof") != data.end()) xdof = data.at("xdof");
    if (data.find("ydof") != data.end()) ydof = data.at("ydof");
    if (data.find("zdof") != data.end()) zdof = data.at("zdof");
    if (data.find("xdoftensor_reduced") != data.end()) xdoftensor_reduced = data.at("xdoftensor_reduced");
    if (data.find("ydoftensor_reduced") != data.end()) ydoftensor_reduced = data.at("ydoftensor_reduced");
    if (data.find("zdoftensor_reduced") != data.end()) zdoftensor_reduced = data.at("zdoftensor_reduced");
    
    // Load grid parameters
    if (data.find("hx") != data.end()) hx = data.at("hx")[0];
    if (data.find("hy") != data.end()) hy = data.at("hy")[0];
    if (data.find("hz") != data.end()) hz = data.at("hz")[0];
    if (data.find("padded_value_size") != data.end()) padded_value_size = static_cast<int>(data.at("padded_value_size")[0]);
    if (data.find("dofs_to_keep") != data.end()) dofs_to_keep = static_cast<uint32_t>(data.at("dofs_to_keep")[0]);
    if (data.find("cells_to_keep") != data.end()) cells_to_keep = static_cast<uint32_t>(data.at("cells_to_keep")[0]);
    if (data.find("local_vals_size") != data.end()) local_vals_size = static_cast<int>(data.at("local_vals_size")[0]);
    
    // Load mappings
    if (data.find("reduced_to_full_map") != data.end()) {
        const auto& vec = data.at("reduced_to_full_map");
        reduced_to_full_map = std::vector<uint32_t>(vec.begin(), vec.end());
    }
    if (data.find("full_to_reduced_map") != data.end()) {
        const auto& vec = data.at("full_to_reduced_map");
        full_to_reduced_map = std::vector<uint32_t>(vec.begin(), vec.end());
    }
    
    // Load skip information
    if (data.find("skip_cell") != data.end()) {
        const auto& vec = data.at("skip_cell");
        skip_cell = std::vector<bool>(vec.begin(), vec.end());
    }
    
    // Load local values
    all_local_vals_map.clear();
    for (const auto& pair : data) {
        if (pair.first.substr(0, 11) == "local_vals_") {
            try {
                int cell_idx = std::stoi(pair.first.substr(11));
                all_local_vals_map[cell_idx] = AlignedPaddedVec(pair.second.begin(), pair.second.end());
            } catch (const std::exception& e) {
                // Skip keys that don't have valid integer indices
                // This can happen if the data contains other keys that start with "local_vals_"
                continue;
            }
        }
    }
}

