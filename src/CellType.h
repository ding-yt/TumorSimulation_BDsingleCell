//
//  CellType.h
//  CancerModel_1
//
//  Created by Yuantong Ding on 11/18/13.
//  Copyright (c) 2013 Yuantong Ding. All rights reserved.
//

#ifndef __CancerModel_1__CellType__
#define __CancerModel_1__CellType__

#include <iostream>
#include <vector>

class Trait
{
protected:
    int name;
    double proliferation_time;
    double death_rate;
    double transition_rate;
    double migration_rate;
    double fittness;
    double point_mutation_rate;
    
};

class Type : public Trait
{
public:
    Type();
    Type (int name_);
    void set_proliferationTime (double time);
    void set_deathRate (double rate);
    void set_mutationRate (double rate);
    void set_migrationRate (double rate);
    void set_fittness (double fit){fittness = fit;};
    void set_point_mutationRate(double rate){point_mutation_rate = rate;};
    
    int get_name() const { return name; }
    double get_proliferation_time() const { return proliferation_time; }
    double get_death_rate() const { return death_rate; }
    double get_mutation_rate() const { return transition_rate; }
    double get_migration_rate() const { return migration_rate; }
    double get_fittness () const {return fittness;}
    double get_point_mutationRate(){return point_mutation_rate;};
    
    
    void show();
    
};


class CellType
{
    std::vector<Type> all_types;
    int type_number;
    
public:
    CellType (int n_type);
    void set_typeNumber(int n_type){type_number = n_type;};
    Type& operator[] (const int nIndex);
    void show();
    double max_P();
    int getTypeNumber(){return type_number;};
    double getMutationRate(int type) {return all_types[type].get_point_mutationRate();};
    void addType(int type_index, double mutationRate, double proliferationRate, double deathRate, double migrationRate, double fitness, double transitionRate){
        type_number ++;
        all_types.push_back(Type(type_index));
        all_types[type_index].set_mutationRate(transitionRate);
        all_types[type_index].set_proliferationTime(proliferationRate);
        all_types[type_index].set_deathRate(deathRate);
        all_types[type_index].set_migrationRate(migrationRate);
        all_types[type_index].set_fittness(fitness);;
        all_types[type_index].set_point_mutationRate(mutationRate);
    }
};


#endif /* defined(__CancerModel_1__CellType__) */
