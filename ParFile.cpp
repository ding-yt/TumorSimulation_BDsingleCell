//
//  ParFile.cpp
//  CancerModel_1
//
//  Created by Yuantong Ding on 11/22/13.
//  Copyright (c) 2013 Yuantong Ding. All rights reserved.
//

#include "ParFile.h"
#include <sstream>

std::map<std::string, double> ParFile::get_parameters (std::string filename)
{
    std::string line;
    parameterFile.open(filename.c_str());
    if (parameterFile.is_open()) {
        while (std::getline(parameterFile, line)) {
            std::istringstream iss(line);
            std::string s, dump;
            double d;
            if (iss >> s >> dump >> d) {
                if (s.find("#") == std::string::npos) {
                    parameters[s] = d;
                }
            }
        }
        parameterFile.close();
    }else{
        std::cout << "Unable to open file "<<filename<<"\n";
    }
    return parameters;
}

std::vector<std::vector<int>> ParFile::get_lineage(std::string filename){
    std::string line;
    std::vector<std::vector<int>> lineages;
    parameterFile.open(filename.c_str());
    if (parameterFile.is_open()) {
        while (std::getline(parameterFile, line)) {
            std::istringstream iss(line);
            int parent;
            std::vector<int> tempLineage;
            while (iss >> parent) {
                tempLineage.push_back(parent);
            }
            lineages.push_back(tempLineage);
        }
        parameterFile.close();
    }else{
        std::cout << "Unable to open file "<<filename<<"\n";
    }
    
    return lineages;
}

void ParFile::show()
{
    for (std::map<std::string, double>::iterator it = parameters.begin(); it != parameters.end(); ++it) {
        std::cout << it->first <<" : "<< it->second <<"\n";
    }
}