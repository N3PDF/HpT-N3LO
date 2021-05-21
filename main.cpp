/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  Main file that computes either the full hadronic or the Mellin
 *                  partonic cross section.
 *
 *        Version:  1.0
 *        Created:  17/02/2021 23:16:58
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tanjona R. Rabemananjara, Roy Stegeman
 *   Organization:  N3PDF
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <exception>

#include "include/CombinedRes.h"
#include "yaml-cpp/yaml.h"
#include "higgs-fo/params.h"
#include "higgs-fo/higgspt.h"
#include "higgs-fo/higgsptpartonic.h"
#include "higgs-fo/partonic.h"
#include "higgs-fo/hadronic.h"
#include "higgs-fo/luminosity.h"


// Exception for wrong inputs
struct err_message : public std::exception {
    const char * what() const throw() {
        return "Wrong Parameters!!";
    }
};


int main(int argc, char* argv[]) {
    LHAPDF::setVerbosity(0);
    YAML::Node node = YAML::LoadFile(argv[1]);

    int inorm = node["inorm"].as<int>();
    int order = node["order"].as<int>();
    int _nf = node["nf"].as<int>();
    int channel = node["channel"].as<int>();
    int scheme = node["scheme"].as<int>();

    double _mh = node["mh"].as<double>();
    double _mur = node["mur"].as<double>();
    double _muf = node["muf"].as<double>();
    double _sroot = node["sroot"].as<double>();
    double ptmin = node["ptmin"].as<double>();
    double ptmax = node["ptmax"].as<double>();
    double ptbin = node["ptbin"].as<double>();
    double nn = node["N"].as<double>();

    /* double y1 = node["y1"].as<double>(); */
    /* double y2 = node["y2"].as<double>(); */

    std::string sectype = node["sectype"].as<std::string>();
    std::string pdfname = node["pdfname"].as<std::string>();
    std::string filename = node["outfile"].as<std::string>();

    std::string ord_fixod[2] = {"_LO", "_NLO"};
    std::string par_chanl[5] = {
        "_gg_channel.dat",
        "_gq_channel.dat",
        "_qq_channel.dat",
        "_qqb_channel.dat",
        "_all_channels.dat"
    };

    try {
        if (order<0 || order>1) throw err_message();
        if (channel<0 || channel>5) throw err_message();
        filename += ord_fixod[order];
        filename += par_chanl[channel];
    } catch(err_message& err) {
        std::cout << err.what() << std::endl;
    }

    // Factors for Born cross-section
    double factor;
    double gf = 1.16637e-5;                         // Fermi Constant
    double gevpb = 3.8937966e8;                     // GeV to pb

    if (inorm == 1) {
        std::cout << "ERROR in inorm!" << std::endl;
        exit(EXIT_FAILURE);                         // TODO: complete implementation!
    } else if (inorm == 0) {
        factor = gf/288./M_PI/sqrt(2.);             // large-top mass limit
    } else {
        std::cout << "ERROR in inorm!" << std::endl;
        exit(EXIT_FAILURE);
    }


    // Initialize PDF and extract alpha_s
    LHAPDF::initPDFSetByName(pdfname);
    double _as = LHAPDF::alphasPDF(_mur);
    double _sigma0 = factor*gevpb*std::pow(_as,2);

    // Define parameters
    PhysParams physparam;
    physparam.nc = 3;
    physparam.nf = _nf;
    physparam.mh = _mh;
    physparam.mur = _mur;
    physparam.muf = _muf;
    physparam.alphas = _as;
    physparam.sroot = _sroot;
    physparam.sigma0 = _sigma0;

    // Init. combined resummation class
    CombinedRes combres(order, channel, pdfname, &physparam);

    // Construct output fie
    std::ofstream output_file(filename);
    output_file << "# Process type         : " << sectype << "\n"
                << "# PDF set name         : " << pdfname << "\n"
                << "# Fixed Order          : " << order   << "\n"
                << "# Partonic channel     : " << channel << "\n"
                << "# Center of M.E. (GeV) : " << _sroot  << "\n"
                << "# Higgs mass (GeV)     : " << _mh     << "\n"
                << "# Renorm. scale (GeV)  : " << _mur    << "\n"
                << "# Fact. scale (GeV)    : " << _muf    << "\n";
    const int space = 16;
    output_file << "# [pt value]" << std::setw(space) << "[dHpt (pb)]"
                << std::setw(space) << "[error (pb)]" << "\n";

    double pt = ptmin;
    std::complex<double> results;
    while (pt <= ptmax) {
        std::complex<double> Ncmpx(nn,0.);
        results = combres.CombinedResExpr(Ncmpx, pt, scheme);

        // Generate some output logs & write to output file
        printf("pt=%e: dHdpt = %e + %e II. \n",
            pt, results.real(), results.imag());
        output_file.setf(std::ios_base::scientific);
        output_file << pt << std::setw(space)
                    << results << "\n";
        output_file.flush();

        pt += ptbin;
    }

    output_file.close();

    return 0;
}
