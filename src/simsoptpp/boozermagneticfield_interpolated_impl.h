#pragma once

#include "boozermagneticfield_interpolated.h"

// Implementation of save/load methods for InterpolatedBoozerField

std::map<std::string, std::map<std::string, std::vector<double>>> InterpolatedBoozerField::get_all_interpolant_data() const {
    std::map<std::string, std::map<std::string, std::vector<double>>> all_data;
    
    // Save data for each interpolant that exists
    if (interp_modB && interp_modB->is_computed()) {
        all_data["modB"] = interp_modB->get_interpolant_data();
    }
    if (interp_dmodBdtheta && interp_dmodBdtheta->is_computed()) {
        all_data["dmodBdtheta"] = interp_dmodBdtheta->get_interpolant_data();
    }
    if (interp_dmodBdzeta && interp_dmodBdzeta->is_computed()) {
        all_data["dmodBdzeta"] = interp_dmodBdzeta->get_interpolant_data();
    }
    if (interp_dmodBds && interp_dmodBds->is_computed()) {
        all_data["dmodBds"] = interp_dmodBds->get_interpolant_data();
    }
    if (interp_G && interp_G->is_computed()) {
        all_data["G"] = interp_G->get_interpolant_data();
    }
    if (interp_I && interp_I->is_computed()) {
        all_data["I"] = interp_I->get_interpolant_data();
    }
    if (interp_iota && interp_iota->is_computed()) {
        all_data["iota"] = interp_iota->get_interpolant_data();
    }
    if (interp_dGds && interp_dGds->is_computed()) {
        all_data["dGds"] = interp_dGds->get_interpolant_data();
    }
    if (interp_dIds && interp_dIds->is_computed()) {
        all_data["dIds"] = interp_dIds->get_interpolant_data();
    }
    if (interp_diotads && interp_diotads->is_computed()) {
        all_data["diotads"] = interp_diotads->get_interpolant_data();
    }
    if (interp_psip && interp_psip->is_computed()) {
        all_data["psip"] = interp_psip->get_interpolant_data();
    }
    if (interp_R && interp_R->is_computed()) {
        all_data["R"] = interp_R->get_interpolant_data();
    }
    if (interp_Z && interp_Z->is_computed()) {
        all_data["Z"] = interp_Z->get_interpolant_data();
    }
    if (interp_nu && interp_nu->is_computed()) {
        all_data["nu"] = interp_nu->get_interpolant_data();
    }
    if (interp_K && interp_K->is_computed()) {
        all_data["K"] = interp_K->get_interpolant_data();
    }
    if (interp_dRdtheta && interp_dRdtheta->is_computed()) {
        all_data["dRdtheta"] = interp_dRdtheta->get_interpolant_data();
    }
    if (interp_dRdzeta && interp_dRdzeta->is_computed()) {
        all_data["dRdzeta"] = interp_dRdzeta->get_interpolant_data();
    }
    if (interp_dRds && interp_dRds->is_computed()) {
        all_data["dRds"] = interp_dRds->get_interpolant_data();
    }
    if (interp_dZdtheta && interp_dZdtheta->is_computed()) {
        all_data["dZdtheta"] = interp_dZdtheta->get_interpolant_data();
    }
    if (interp_dZdzeta && interp_dZdzeta->is_computed()) {
        all_data["dZdzeta"] = interp_dZdzeta->get_interpolant_data();
    }
    if (interp_dZds && interp_dZds->is_computed()) {
        all_data["dZds"] = interp_dZds->get_interpolant_data();
    }
    if (interp_dnudtheta && interp_dnudtheta->is_computed()) {
        all_data["dnudtheta"] = interp_dnudtheta->get_interpolant_data();
    }
    if (interp_dnudzeta && interp_dnudzeta->is_computed()) {
        all_data["dnudzeta"] = interp_dnudzeta->get_interpolant_data();
    }
    if (interp_dnuds && interp_dnuds->is_computed()) {
        all_data["dnuds"] = interp_dnuds->get_interpolant_data();
    }
    if (interp_dKdtheta && interp_dKdtheta->is_computed()) {
        all_data["dKdtheta"] = interp_dKdtheta->get_interpolant_data();
    }
    if (interp_dKdzeta && interp_dKdzeta->is_computed()) {
        all_data["dKdzeta"] = interp_dKdzeta->get_interpolant_data();
    }
    if (interp_K_derivs && interp_K_derivs->is_computed()) {
        all_data["K_derivs"] = interp_K_derivs->get_interpolant_data();
    }
    if (interp_nu_derivs && interp_nu_derivs->is_computed()) {
        all_data["nu_derivs"] = interp_nu_derivs->get_interpolant_data();
    }
    if (interp_R_derivs && interp_R_derivs->is_computed()) {
        all_data["R_derivs"] = interp_R_derivs->get_interpolant_data();
    }
    if (interp_Z_derivs && interp_Z_derivs->is_computed()) {
        all_data["Z_derivs"] = interp_Z_derivs->get_interpolant_data();
    }
    if (interp_modB_derivs && interp_modB_derivs->is_computed()) {
        all_data["modB_derivs"] = interp_modB_derivs->get_interpolant_data();
    }
    
    return all_data;
}

void InterpolatedBoozerField::set_all_interpolant_data(const std::map<std::string, std::map<std::string, std::vector<double>>>& data) {
    // Load data for each interpolant
    for (const auto& pair : data) {
        const std::string& quantity = pair.first;
        const std::map<std::string, std::vector<double>>& interpolant_data = pair.second;
        
        if (quantity == "modB" && interp_modB) {
            interp_modB->set_interpolant_data(interpolant_data);
        } else if (quantity == "dmodBdtheta" && interp_dmodBdtheta) {
            interp_dmodBdtheta->set_interpolant_data(interpolant_data);
        } else if (quantity == "dmodBdzeta" && interp_dmodBdzeta) {
            interp_dmodBdzeta->set_interpolant_data(interpolant_data);
        } else if (quantity == "dmodBds" && interp_dmodBds) {
            interp_dmodBds->set_interpolant_data(interpolant_data);
        } else if (quantity == "G" && interp_G) {
            interp_G->set_interpolant_data(interpolant_data);
        } else if (quantity == "I" && interp_I) {
            interp_I->set_interpolant_data(interpolant_data);
        } else if (quantity == "iota" && interp_iota) {
            interp_iota->set_interpolant_data(interpolant_data);
        } else if (quantity == "dGds" && interp_dGds) {
            interp_dGds->set_interpolant_data(interpolant_data);
        } else if (quantity == "dIds" && interp_dIds) {
            interp_dIds->set_interpolant_data(interpolant_data);
        } else if (quantity == "diotads" && interp_diotads) {
            interp_diotads->set_interpolant_data(interpolant_data);
        } else if (quantity == "psip" && interp_psip) {
            interp_psip->set_interpolant_data(interpolant_data);
        } else if (quantity == "R" && interp_R) {
            interp_R->set_interpolant_data(interpolant_data);
        } else if (quantity == "Z" && interp_Z) {
            interp_Z->set_interpolant_data(interpolant_data);
        } else if (quantity == "nu" && interp_nu) {
            interp_nu->set_interpolant_data(interpolant_data);
        } else if (quantity == "K" && interp_K) {
            interp_K->set_interpolant_data(interpolant_data);
        } else if (quantity == "dRdtheta" && interp_dRdtheta) {
            interp_dRdtheta->set_interpolant_data(interpolant_data);
        } else if (quantity == "dRdzeta" && interp_dRdzeta) {
            interp_dRdzeta->set_interpolant_data(interpolant_data);
        } else if (quantity == "dRds" && interp_dRds) {
            interp_dRds->set_interpolant_data(interpolant_data);
        } else if (quantity == "dZdtheta" && interp_dZdtheta) {
            interp_dZdtheta->set_interpolant_data(interpolant_data);
        } else if (quantity == "dZdzeta" && interp_dZdzeta) {
            interp_dZdzeta->set_interpolant_data(interpolant_data);
        } else if (quantity == "dZds" && interp_dZds) {
            interp_dZds->set_interpolant_data(interpolant_data);
        } else if (quantity == "dnudtheta" && interp_dnudtheta) {
            interp_dnudtheta->set_interpolant_data(interpolant_data);
        } else if (quantity == "dnudzeta" && interp_dnudzeta) {
            interp_dnudzeta->set_interpolant_data(interpolant_data);
        } else if (quantity == "dnuds" && interp_dnuds) {
            interp_dnuds->set_interpolant_data(interpolant_data);
        } else if (quantity == "dKdtheta" && interp_dKdtheta) {
            interp_dKdtheta->set_interpolant_data(interpolant_data);
        } else if (quantity == "dKdzeta" && interp_dKdzeta) {
            interp_dKdzeta->set_interpolant_data(interpolant_data);
        } else if (quantity == "K_derivs" && interp_K_derivs) {
            interp_K_derivs->set_interpolant_data(interpolant_data);
        } else if (quantity == "nu_derivs" && interp_nu_derivs) {
            interp_nu_derivs->set_interpolant_data(interpolant_data);
        } else if (quantity == "R_derivs" && interp_R_derivs) {
            interp_R_derivs->set_interpolant_data(interpolant_data);
        } else if (quantity == "Z_derivs" && interp_Z_derivs) {
            interp_Z_derivs->set_interpolant_data(interpolant_data);
        } else if (quantity == "modB_derivs" && interp_modB_derivs) {
            interp_modB_derivs->set_interpolant_data(interpolant_data);
        }
    }
}

std::map<std::string, bool> InterpolatedBoozerField::get_status_flags() const {
    std::map<std::string, bool> flags;
    flags["status_modB"] = status_modB;
    flags["status_dmodBdtheta"] = status_dmodBdtheta;
    flags["status_dmodBdzeta"] = status_dmodBdzeta;
    flags["status_dmodBds"] = status_dmodBds;
    flags["status_G"] = status_G;
    flags["status_I"] = status_I;
    flags["status_iota"] = status_iota;
    flags["status_dGds"] = status_dGds;
    flags["status_dIds"] = status_dIds;
    flags["status_diotads"] = status_diotads;
    flags["status_psip"] = status_psip;
    flags["status_R"] = status_R;
    flags["status_Z"] = status_Z;
    flags["status_nu"] = status_nu;
    flags["status_K"] = status_K;
    flags["status_dRdtheta"] = status_dRdtheta;
    flags["status_dRdzeta"] = status_dRdzeta;
    flags["status_dRds"] = status_dRds;
    flags["status_dZdtheta"] = status_dZdtheta;
    flags["status_dZdzeta"] = status_dZdzeta;
    flags["status_dZds"] = status_dZds;
    flags["status_dnudtheta"] = status_dnudtheta;
    flags["status_dnudzeta"] = status_dnudzeta;
    flags["status_dnuds"] = status_dnuds;
    flags["status_dKdtheta"] = status_dKdtheta;
    flags["status_dKdzeta"] = status_dKdzeta;
    flags["status_K_derivs"] = status_K_derivs;
    flags["status_R_derivs"] = status_R_derivs;
    flags["status_Z_derivs"] = status_Z_derivs;
    flags["status_nu_derivs"] = status_nu_derivs;
    flags["status_modB_derivs"] = status_modB_derivs;
    return flags;
}

void InterpolatedBoozerField::set_status_flags(const std::map<std::string, bool>& flags) {
    if (flags.find("status_modB") != flags.end()) status_modB = flags.at("status_modB");
    if (flags.find("status_dmodBdtheta") != flags.end()) status_dmodBdtheta = flags.at("status_dmodBdtheta");
    if (flags.find("status_dmodBdzeta") != flags.end()) status_dmodBdzeta = flags.at("status_dmodBdzeta");
    if (flags.find("status_dmodBds") != flags.end()) status_dmodBds = flags.at("status_dmodBds");
    if (flags.find("status_G") != flags.end()) status_G = flags.at("status_G");
    if (flags.find("status_I") != flags.end()) status_I = flags.at("status_I");
    if (flags.find("status_iota") != flags.end()) status_iota = flags.at("status_iota");
    if (flags.find("status_dGds") != flags.end()) status_dGds = flags.at("status_dGds");
    if (flags.find("status_dIds") != flags.end()) status_dIds = flags.at("status_dIds");
    if (flags.find("status_diotads") != flags.end()) status_diotads = flags.at("status_diotads");
    if (flags.find("status_psip") != flags.end()) status_psip = flags.at("status_psip");
    if (flags.find("status_R") != flags.end()) status_R = flags.at("status_R");
    if (flags.find("status_Z") != flags.end()) status_Z = flags.at("status_Z");
    if (flags.find("status_nu") != flags.end()) status_nu = flags.at("status_nu");
    if (flags.find("status_K") != flags.end()) status_K = flags.at("status_K");
    if (flags.find("status_dRdtheta") != flags.end()) status_dRdtheta = flags.at("status_dRdtheta");
    if (flags.find("status_dRdzeta") != flags.end()) status_dRdzeta = flags.at("status_dRdzeta");
    if (flags.find("status_dRds") != flags.end()) status_dRds = flags.at("status_dRds");
    if (flags.find("status_dZdtheta") != flags.end()) status_dZdtheta = flags.at("status_dZdtheta");
    if (flags.find("status_dZdzeta") != flags.end()) status_dZdzeta = flags.at("status_dZdzeta");
    if (flags.find("status_dZds") != flags.end()) status_dZds = flags.at("status_dZds");
    if (flags.find("status_dnudtheta") != flags.end()) status_dnudtheta = flags.at("status_dnudtheta");
    if (flags.find("status_dnudzeta") != flags.end()) status_dnudzeta = flags.at("status_dnudzeta");
    if (flags.find("status_dnuds") != flags.end()) status_dnuds = flags.at("status_dnuds");
    if (flags.find("status_dKdtheta") != flags.end()) status_dKdtheta = flags.at("status_dKdtheta");
    if (flags.find("status_dKdzeta") != flags.end()) status_dKdzeta = flags.at("status_dKdzeta");
    if (flags.find("status_K_derivs") != flags.end()) status_K_derivs = flags.at("status_K_derivs");
    if (flags.find("status_R_derivs") != flags.end()) status_R_derivs = flags.at("status_R_derivs");
    if (flags.find("status_Z_derivs") != flags.end()) status_Z_derivs = flags.at("status_Z_derivs");
    if (flags.find("status_nu_derivs") != flags.end()) status_nu_derivs = flags.at("status_nu_derivs");
    if (flags.find("status_modB_derivs") != flags.end()) status_modB_derivs = flags.at("status_modB_derivs");
}
