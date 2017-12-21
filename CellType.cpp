//
//  CellType.cpp
//  CancerModel_1
//
//  Created by Yuantong Ding on 11/18/13.
//  Copyright (c) 2013 Yuantong Ding. All rights reserved.
//

#include "CellType.h"
//#include <iostream>

Type::Type (int name_)
{
    name = name_;
    transition_rate = 0;
    migration_rate = 0;
    proliferation_time = 0;
    death_rate = 0;
    point_mutation_rate = 0;
    
}

void Type::set_proliferationTime (double time)
{
    proliferation_time = time;
}
void Type::set_deathRate (double rate)
{
    death_rate = rate;
}
void Type::set_mutationRate (double rate)
{
    transition_rate = rate;
}
void Type::set_migrationRate (double rate)
{
    migration_rate = rate;
}


void Type::show()
{
    std::cout <<"Type name: "<<name<<"\n";
    std::cout <<"Cell prolifertation: "<<proliferation_time<<"\n";
    std::cout <<"Cell death rate: "<<death_rate<<"\n";
    std::cout <<"Cell migration rate: "<<migration_rate<<"\n";
    std::cout <<"Cell transition rate: "<<transition_rate<<"\n";
    std::cout <<"Cell fittness:"<<fittness<<"\n";
    std::cout <<"Cell point mutation rate:"<<point_mutation_rate<<"\n";
}


CellType::CellType (int n_type)
{
    type_number = n_type;
    for (int i=0; i<n_type; i++) {
        all_types.push_back(Type(i));
    }
}

Type& CellType::operator[] (const int nIndex)
{
    return all_types[nIndex];
}

void CellType::show()
{
    for (int i=0 ;i < type_number;  i++) {
        all_types[i].show();
        std::cout <<"\n";
    }
}

double CellType::max_P(){
    double p=0;
    for (int i=0;i<all_types.size();i++){
        if (p<all_types[i].get_proliferation_time()) {
            p=all_types[i].get_proliferation_time();
        }
    }
    return p;
}