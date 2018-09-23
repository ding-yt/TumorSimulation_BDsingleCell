//
//  GenealogyNode.h
//  CancerModel_cleanup
//
//  Created by Yuantong Ding on 2/26/14.
//  Copyright (c) 2014 Yuantong Ding. All rights reserved.
//

#ifndef __CancerModel_cleanup__GenealogyNode__
#define __CancerModel_cleanup__GenealogyNode__

#include <iostream>
#include <algorithm>
#include <cassert>
#include <vector>
#include <iterator>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <cstring>
#include <stack>
#include <list>
#include "CellType.h"

///////////////////////////////////////////////////////////////////////////////
// GenealogyNode
/**
 * A single node of a genealogical tree.
 *
 * Represents a single node of a genealogy of a neutral marker, with pointers
 * to its parent node. Uses reference counting to tracking references to self,
 * and deletes self when no other node points to self as parent.
 */

//typedef unsigned int           CellIndexType;
typedef int           CellIndexType;

class GenealogyNode {
    
private:
    
    /**
     * Pointer to the parent of this node (<code>NULL</code> if no
     * children).
     */
    GenealogyNode *     parent_;
    
    /**
     * Number of objects that point to or reference this object (including
     * this object itself).
     */
    unsigned            reference_count_;
    
    /**
     * Position of containing organism: for reconstruction.
     */
    CellIndexType       cell_index_;
    
    /**
     * stage of node, alive or dead.
     */
    bool alive_;
    
    /**
     * Birth or death time of the node.
     */
    double birth_time_;
    double death_time_;
    
    /**
     * cell type at this node.
     */
    int type_;
    
//    std::pair<int, int> _loc;
    
    /**
     * static cell type information.
     */
    static CellType celltypes_;
    
    static int _count;
    
    
public:
    GenealogyNode();
    
    GenealogyNode(int type, double birth_time);
    
    ~GenealogyNode();
    
    void clear();
    
    unsigned reference_count() const;

    GenealogyNode * get_parent();
    
    void set_parent(GenealogyNode * parent);
    
    CellIndexType get_cell_index();
    
    void set_cell_index(CellIndexType cell_index);
    
    void die(double death_time);
    
    bool isAlive();
    
    void set_birth_time(double time);
    
    void set_death_time(double time);
    
    double get_birth_time();
    
    double get_death_time();
    
    void set_type(int type);
    
    int get_type();
    
    void set_cell_type(int type){type_ = type;};
    
    double get_proliferation_time() const { return celltypes_[type_].get_proliferation_time(); }
    
    double get_death_rate() const { return celltypes_[type_].get_death_rate(); }
    
    double get_mutation_rate() const { return celltypes_[type_].get_mutation_rate(); }
    
    double get_migration_rate() const { return celltypes_[type_].get_migration_rate(); }
    
    double get_birthTime() const { return birth_time_; }
    
    double get_deathTime() const { return death_time_; }
    
    double get_fittness() const { return celltypes_[type_].get_fittness(); }
    
    double get_maxP() { return celltypes_.max_P(); }
    
    static void show_allCellType(){celltypes_.show();};
    
    static void set_counter(int i){_count = i;};
    
    static void set_allCellType(std::map<std::string, double> parameters){
        
        int cell_type_number = parameters["cell_types"];
        CellType allCellType(cell_type_number);
        double max_p=0;
        for (int i=0; i<cell_type_number; i++) {
            std::string prefix = "type_";
            std::ostringstream s;
            s << i;
            std::string proliferation = prefix + s.str() + "_proliferation";
            if (max_p < parameters[proliferation]){
                max_p = parameters[proliferation];
            }
        }

        for (int i=0; i<cell_type_number; i++) {
            std::string prefix = "type_";
            std::ostringstream s;
            s << i;
            std::string proliferation = prefix + s.str() + "_proliferation";
            std::string mutation = prefix + s.str() + "_transition_rate";
            std::string migration = prefix + s.str() + "_migration_rate";
            std::string death = prefix + s.str() + "_death_rate";
            std::string fitness = prefix + s.str() + "_fitness";
            std::string point_mutation_rate = prefix + s.str() + "_point_mutation_rate";
            allCellType[i].set_proliferationTime(parameters[proliferation]);
            allCellType[i].set_mutationRate(parameters[mutation]);
            allCellType[i].set_migrationRate(parameters[migration]);
            allCellType[i].set_deathRate(parameters[death]);
            allCellType[i].set_fittness(parameters[fitness]);
            allCellType[i].set_point_mutationRate(parameters[point_mutation_rate]);
        }
        celltypes_ = allCellType;
        
    };
    
    void remove();
    
    static int cellTypeCount() {return celltypes_.getTypeNumber();};
    
    static double getMutationRate(int type) { return celltypes_[type].get_point_mutationRate();};
    
    static double getProliferationRate(int type) { return celltypes_[type].get_proliferation_time();};
    
    static double getDeathRate(int type) { return celltypes_[type].get_death_rate();};
    
    static double getMigrationRate(int type) { return celltypes_[type].get_migration_rate();};
    
    static double getFittness(int type) { return celltypes_[type].get_fittness();};
    
    static double getTransitionRate(int type) { return celltypes_[type].get_mutation_rate();};
    
    static bool existType(int type) {
        if (celltypes_.getTypeNumber()> type) {
            return true;
        } else {
            return false;
        }
    };
    
    static void addType(int type, double mutationRate, double proliferationRate, double deathRate, double migrationRate, double fitness, double transitionRate){
        celltypes_.addType(type, mutationRate, proliferationRate, deathRate, migrationRate, fitness, transitionRate);
    }
    
private:
    
    /** Copy constructor (disabled by private scoping). */
    GenealogyNode(const GenealogyNode& );
    
    /** Assignment constructor (disabled by private scoping). */
    const GenealogyNode& operator=(const GenealogyNode&);
    
    void decrement_count();
    
    void increment_count();
    
};


// GenealogyNode
///////////////////////////////////////////////////////////////////////////////

#endif /* defined(__CancerModel_cleanup__GenealogyNode__) */
