//
//  GenealogyNode.cpp
//  CancerModel_cleanup
//
//  Created by Yuantong Ding on 2/26/14.
//  Copyright (c) 2014 Yuantong Ding. All rights reserved.
//

#include "GenealogyNode.h"
CellType GenealogyNode::celltypes_ =CellType(1);
int GenealogyNode::_count = 0;

using namespace std;
/** Constructs a node with no antecedents. */
GenealogyNode::GenealogyNode(){
    parent_ = NULL;
    reference_count_ = 1;
    cell_index_ = _count;
    _count++;
    alive_=true;
    birth_time_=0;
    death_time_=-1;
    type_=0;

 }


GenealogyNode::GenealogyNode(int type, double birth_time){
    this->parent_= NULL;
    this->reference_count_ = 1;
    this->cell_index_ = _count;
    this->_count++;
    this->alive_ = true;
    this->birth_time_ = birth_time_;
    this->death_time_ = -1;
    this->type_ = type;

}

/**
 * Ensures all pointers from and to this object are nulled out before
 * winking out of existence.
 */
GenealogyNode::~GenealogyNode() {
    //std::cout<<"deleting: "<<cell_index_<<"\n";
    this->clear();
}

void GenealogyNode::clear() {
    assert(this->parent_ != this);
    assert(this->reference_count() == 1);
    //cout<<"clear\n";
    if (this->parent_) {
        this->parent_->decrement_count();
        //this->parent_ = NULL;
    }
}


/**
 * Registers one less reference to this node.
 *
 * Decrements the count of references to this node (number of objects
 * that point to this node). If this goes to 1 (which means that the
 * only object that references this node object is this node object
 * itself), then this deletes itself.
 */
void GenealogyNode::decrement_count() {
    if (this->reference_count_ == 1) {
        delete this;
    } else {
        this->reference_count_ -= 1;
    
    if (this->reference_count_ == 1) {
        delete this;
    }
    }
}

/**
 * Registers one more reference to this node.
 *
 * Notes that one more object points to or otherwise references this
 * node object.
 */
void GenealogyNode::increment_count() {
    this->reference_count_ += 1;
}

/**
 * Set a node to death
 * If the node is a leaf (_count=1), delete it
 *
 */
void GenealogyNode::die(double death_time){
    alive_=false;
    death_time_ = death_time;
    if (this->reference_count()==1) {
        //cout<<"delete\n";
        delete this;
    }
}

bool GenealogyNode::isAlive(){
    return alive_;
}

void GenealogyNode::set_birth_time(double time){
    birth_time_ = time;
};

void GenealogyNode::set_death_time(double time){
    death_time_ = time;
    alive_ = false; //modify: added 9/8/15
};

double GenealogyNode::get_birth_time(){
    return birth_time_;
};

double GenealogyNode::get_death_time(){
    return death_time_;
};

void GenealogyNode::set_type( int type){
    type_ = type;
};

int GenealogyNode::get_type(){
    return type_;
};

/**
 * Returns the cell index of this node.
 */
CellIndexType GenealogyNode::get_cell_index() {
    return this->cell_index_;
}

////////////////////???????????????????????????
/**
 * Sets the cell index.
 * @param cell_index    index of cell occupied by current organism.
 */
void GenealogyNode::set_cell_index(CellIndexType cell_index) {
    this->cell_index_ = cell_index;
}

/**
 * Returns pointer to parent node.
 *
 * @return      pointer to parent node
 */
GenealogyNode * GenealogyNode::get_parent() {
    return this->parent_;
}

/**
 * Sets pointer to parent node.
 *
 * @param parent    pointer to parent node
 */
void GenealogyNode::set_parent(GenealogyNode * parent) {
    assert(parent != this);
    if (this->parent_ != NULL){
        this->parent_->decrement_count();
    }
    this->parent_ = parent;
    
    if (this->parent_ != NULL) {
        this->parent_->increment_count();
    }
}


///////////////?????
/**
 * Returns reference count.
 */
unsigned GenealogyNode::reference_count() const {
    return this->reference_count_;
}




