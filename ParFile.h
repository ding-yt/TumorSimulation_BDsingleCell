//
//  ParFile.h
//  CancerModel_1
//
//  Created by Yuantong Ding on 11/22/13.
//  Copyright (c) 2013 Yuantong Ding. All rights reserved.
//

#ifndef __CancerModel_1__ParFile__
#define __CancerModel_1__ParFile__

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

class ParFile
{
    std::map<std::string, double> parameters;
    std::ifstream parameterFile;
    std::string filename;
    
public:
    std::map<std::string, double> get_parameters (std::string filename);
    std::vector<std::vector<int>> get_lineage(std::string filename);
    void show();
    
};


#endif /* defined(__CancerModel_1__ParFile__) */
