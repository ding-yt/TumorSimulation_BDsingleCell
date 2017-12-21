//
//  Topology.cpp
//  NodeTester
//
//  Created by Yuantong Ding on 1/27/15.
//  Copyright (c) 2015 Yuantong Ding. All rights reserved.
//

#include "Topology.h"
#include <cassert>
#include <math.h>

Nodetree::Nodetree()
{
    _parent = NULL;
    _index = -1;
    _type = -1;
    _weight = 1;
    _name = "u";
    _sequence = "";
    _locX = -1;
    _locY = -1;
    
}

Nodetree::Nodetree(int index, int type, Nodetree * parent)
{
    _parent = parent;
    _index = index;
    _type = type;
    _weight = 1;
    _name = "u";
    _sequence = "";
    _locX = -1;
    _locY = -1;
}

Nodetree::~Nodetree(){
    
    for (int i=0; i<_children.size(); i++) {
        delete _children[i];
    }
    
}

void Nodetree::addName(int type, int x, int y){
    _name = "t";
    _name.append(std::to_string(static_cast<long long>(type)));
    _name.append("_");
    _name.append(std::to_string(static_cast<long long>(x)));
    _name.append("_");
    _name.append(std::to_string(static_cast<long long>(y)));
}

Topology::Topology():
_rng(ginkgo::RandomNumberGenerator::get_instance())
{
    _root = NULL;
    
}

Topology::~Topology(){
    destoryTree(_root);
}

void Topology::destoryTree(Nodetree * node){
    if (node != NULL){
        for (int i=0; i<node->_children.size(); i++) {
            delete node->_children[i];
        }
        delete node;
    }
}

Nodetree * Topology::searchNode(int index){
    //std::cout<<"root index:"<<_root->_index<<"\n";
    return searchNode(index, _root);
}

Nodetree * Topology::searchNode(int index, Nodetree * startNode){
    //std::cout<<"start Node index :"<<startNode->_index<<"\n";
    if (startNode==NULL) {
        return NULL;
    }
    if (startNode->_index==index){
        //std::cout<<"found:"<<index<<"\n";
        //std::cout<<startNode->_sequence<<"\n";
        return startNode;
    }else{
        if (startNode->_children.size()>0) {
            for (int i=0; i<startNode->_children.size(); i++) {
                Nodetree * n = searchNode(index,startNode->_children[i]);
                if (n != NULL) {
                    return n;
                }
            }
        }
        return NULL;
    }
}

Nodetree * Topology::searchChildren(int index, Nodetree * parentNode){
    if (parentNode == NULL) {
        return NULL;
    }
    for (int i=0; i<parentNode->_children.size(); i++) {
        if (parentNode->_children[i]->_index == index) {
            return parentNode->_children[i];
        }
    }
    return NULL;
    
}


Nodetree * Topology::addNode(Nodetree * parentNode, int index, int type){
    if (parentNode == NULL) {
        parentNode = new Nodetree(index, type, NULL);
        return parentNode;
    }
    Nodetree * children = searchChildren(index, parentNode);
    if (children == NULL) {
        Nodetree * newchild = new Nodetree(index,type,parentNode);
        parentNode->_children.push_back(newchild);
        return parentNode->_children[parentNode->_children.size()-1];
    }else{
        return children;
    }
    
}

void Topology::addLineage(std::vector<CellIndexType> lineage){
    //std::cout <<lineage.size()<<"\n";
    int locX = lineage[lineage.size()-2];
    int locY = lineage[lineage.size()-1];;
    lineage.pop_back();
    lineage.pop_back();
    if (lineage.size()<2) {
        return;
    }
    int endIndex = lineage.size()-1;
    //Nodetree * current = addNode(_root, lineage[endIndex-1], lineage[endIndex]);
    //_root = current;
    Nodetree * current = searchNode(lineage[endIndex-1]);
    for (int i=endIndex-2; i>=0; i-=2) {
        //std::cout<<"add cell:"<<lineage[i-1]<<" type: "<<lineage[i]<<"\n";
        current = addNode(current, lineage[i-1], lineage[i]);
        //        if (i==1) {
        //            current->addName(lineage[i], lineage[endIndex-2], lineage[endIndex-1]);
        //        }
    }
    current->_locX = locX;
    current->_locY = locY;
}

void Topology::printNWK(){
    print(_root);
    std::cout<<"\n";
}

void Topology::print(Nodetree * node){
    if (node ==NULL) {
        std::cout<<"empty tree\n";
        return;
    }
    if (node->_children.size()==0) {
        
        std::cout<<"t"<<node->_type<<"_"<<node->_index<<":"<<node->_weight;
        //std::cout<<node->_index<<":"<<node->_type<<":"<<node->_weight;
        return;
    }
    if (node->_children.size()==1) {
        print(node->_children[0]);
//        std::cout<<"t"<<node->_type<<"_"<<node->_index<<":"<<node->_weight;
//        //std::cout<<node->_index<<":"<<node->_type<<":"<<node->_weight;
//        return;
    }else{
    std::cout<<"(";
    for (int i=0; i<node->_children.size(); i++) {
        print(node->_children[i]);
        if (i!=node->_children.size()-1) {
            std::cout<<",";
        }
        
        //        }
    }
    //std::cout<<")"<<node->_index<<":"<<node->_type<<":"<<node->_weight;
    std::cout<<")"<<node->_type<<":"<<node->_index<<":"<<node->_weight;
    }
}

void Topology::addRoot(int index, int type){
    _root = new Nodetree(index,type,NULL);
}

void Topology::addChild(Nodetree *child){
    if (_root==NULL) {
        std::cout<<"no root!\n";
        return;
    }else{
        _root->_children.push_back(child);
    }
}

void Topology::compress(){
    compress(_root);
}

void Topology::compress_noSingleton(){
    compress_noSingleton(_root);
}

void Topology::compress(Nodetree * node){
    if (node->_children.size()==0) {
        return;
    }
    if (node->_children.size()==1 && node->_parent != NULL && node->_type == node->_children[0]->_type) {
        Nodetree * parent = node->_parent;
        Nodetree * child = node->_children[0];
        child->_weight += node->_weight;
        child->_parent = parent;
        for (int i=0; i<parent->_children.size(); i++) {
            if (parent->_children[i]->_index == node->_index) {
                parent->_children[i] = child;
            }
        }
        node->_parent = NULL;
        node->_children[0] = NULL;
        delete node;
        compress(child);
    }else{
        for (int i=0; i<node->_children.size(); i++) {
            compress(node->_children[i]);
        }
    }
}

void Topology::compress_noSingleton(Nodetree * node){
    if (node->_children.size()==0) {
        return;
    }
    if (node->_children.size()==1 && node->_parent != NULL ) {
        Nodetree * parent = node->_parent;
        Nodetree * child = node->_children[0];
        child->_weight += node->_weight;
        child->_parent = parent;
        for (int i=0; i<parent->_children.size(); i++) {
            if (parent->_children[i]->_index == node->_index) {
                parent->_children[i] = child;
            }
        }
        node->_parent = NULL;
        node->_children[0] = NULL;
        delete node;
        compress(child);
    }else{
        for (int i=0; i<node->_children.size(); i++) {
            compress(node->_children[i]);
        }
    }
}


void Topology::simulateSeq(std::vector<double> rate,int length){
    _root->_sequence = randomBase(length);
    int total_mutations = 0;
    mutateSeq(_root, rate, &total_mutations);
    if (length>total_mutations) {
        chopSeq(_root,total_mutations);
    }
    std::cout<<"total mutations:"<<total_mutations<<"\n";
}

void Topology::chopSeq(Nodetree * node,int length){
    if (node->_children.size()==0) {
        //std::cout<<"seq "<<node->_index<<" previous:"<<node->_sequence.length();
        node->_sequence = node->_sequence.substr(0,length);
        //std::cout<<" after chop: "<<node->_sequence.length()<<"\n";
        return;
    }else{
        //std::cout<<"seq "<<node->_index<<" previous:"<<node->_sequence.length();
        node->_sequence = node->_sequence.substr(0,length);
        //std::cout<<" after chop: "<<node->_sequence.length()<<"\n";
        for (int i=0; i<node->_children.size(); i++) {
            chopSeq(node->_children[i], length);
        }
    }
    
}

void Topology::mutateSeq(Nodetree * node, std::vector<double> rates, int * count){
    //std::cout <<">"<<node->_index<<"_t"<<node->_type<<"_w"<<node->_weight<<"\n";

    if (node->_children.size()==0) {
        //std::cout <<">"<<node->_index<<"_t"<<node->_type<<"_w"<<node->_weight<<"\n"<<node->_sequence<<"\n";
        
        //return count;
        return;
    }
    if (*count>=node->_sequence.size()) {
        std::cout<<"warning: initial seq length too short!\n";
        simulateSeq(rates, node->_sequence.size()*2);
        return;
    }else{
    
    for (int i=0; i<node->_children.size(); i++) {
        node->_children[i]->_sequence = node->_sequence;
        //std::cout <<count<<"\n";
        float rate;
        if (node->_children[i]->_type<0) {
            rate = node->_children[i]->_weight * rates[0];
        }else{
            rate = node->_children[i]->_weight * rates[node->_children[i]->_type];
        }
        //int m = _rng.poisson(rate);
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        int m = round(rate);
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //std::cout<<"p"<<node->_index<<"\tc"<<node->_children[i]->_index<<"\tm"<<m<<"\n";
        if (*count+m>node->_sequence.size()) {
            std::cout<<node->_sequence.size()<<" warning: initial seq length too short!\n";
            simulateSeq(rates, node->_sequence.size()*2);
            return;
        }
        for (int j=0; j<m; j++) {
            //std::cout <<node->_children[i]->_index<<"\n";
            node->_children[i]->_sequence[*count+j-1] = pointMutate(node->_children[i]->_sequence[*count+j]);
        }
        *count += m;
    }
    for (int i=0; i<node->_children.size(); i++) {
        mutateSeq(node->_children[i], rates, count);
    }
    }
};

std::string Topology::randomBase(int n){
    long r;
    std::vector<std::string> bases{"A","T","C","G"};
    std::string s;
    for (int i=0;i<n; i++) {
        r = _rng.uniform_int(0, 3);
        s += bases[r];
    }
    return s;
}
char Topology::pointMutate(char s){
    std::vector<char> bases;
    if (s == 'A'){
        bases = {'T','C','G'};
    }else if (s=='T'){
        bases = {'A','C','G'};
    }else if (s=='C'){
        bases = {'A','T','G'};
    }else if (s=='G'){
        bases = {'A','C','T'};
    }else{
        std::cout<<"warning: unknown base "<<s<<"\n";
    }
    long r = _rng.uniform_int(0, 2);
    return bases[r];
    
}

std::vector<std::string> Topology::printSeq(){
    std::vector<std::string> seqs;
    seqs = getSeq(_root, seqs);
    
}

std::vector<std::string> Topology::getSeq(Nodetree *node, std::vector<std::string> seq){
    
    if (node->_children.size()==0) {
        std::string s = std::to_string(static_cast<long long>(node->_index)) +"_t"+std::to_string(static_cast<long long>(node->_type))+":"+node->_sequence;
        seq.push_back(s);
        return seq;
    }
    
    for (int i=0; i<node->_children.size(); i++) {
        getSeq(node->_children[i],seq);
    }
}

std::string Topology::getSeq(Nodetree *node, int index){
    if (node==NULL) {
        std::string s = "";
        return s;
    }
    if (node->_index==index) {
        std::cout<<"find "<<index<<"\n";
        std::cout<<node->_sequence;
        return node->_sequence;
    }
    if (node->_children.size()>0) {
        for (int i=0; i<node->_children.size();i++){
            std::string n = getSeq(node->_children[i],index);
//            if (n!="") {
//                return n;
//            }
        }
        return "";
    }
    //return "";
}

std::string Topology::getSequence(int index){
        Nodetree * n = searchNode(index, _root);
        if (n==NULL) {
            return "";
        }else{
            return n->_sequence;
        }
    //return getSeq(_root, index);
}

std::string Topology::getNWK(){
    std::string s = "";
    std::string * lineage = &s;
    lineage = getNWK(_root, lineage);
    return *lineage;
}

std::string Topology::getNWK_noSingle(){
    std::cout<<"using getNWK_noSingle here\n\n";
    std::string s = "";
    std::string * lineage = &s;
    lineage = getNWK_noSingle(_root, lineage);
    return *lineage;
}

std::string * Topology::getNWK(Nodetree * node, std::string *lineage){
    if (node ==NULL) {
        return lineage;
    }
    if (node->_children.size()==0) {
        //std::string temp = "t"+std::to_string(node->_type)+"_"+std::to_string(node->_index)+":"+std::to_string(node->_weight);
        std::string temp = "t"+std::to_string(static_cast<long long>(node->_type))+"_"+std::to_string(static_cast<long long>(node->_locX))+"_"+std::to_string(static_cast<long long>(node->_locY))+":"+std::to_string(static_cast<long long>(node->_weight));
        *lineage += temp;
        return lineage;
    }
    //std::cout<<"(";
    *lineage += "(";
    for (int i=0; i<node->_children.size(); i++) {
        //print(node->_children[i]);
        lineage = getNWK(node->_children[i],lineage);
        if (i!=node->_children.size()-1) {
            //std::cout<<",";
            *lineage += ",";
        }
    }
    //std::cout<<")"<<node->_index<<":"<<node->_weight;
    *lineage += ")"+std::to_string(static_cast<long long>(node->_index))+":"+std::to_string(static_cast<long long>(node->_weight));
    return lineage;
}

std::string * Topology::getNWK_noSingle(Nodetree * node, std::string *lineage){
    std::cout<<"geting non single\n";
    if (node ==NULL) {
        return lineage;
    }
    if (node->_children.size()==0) {
        //std::string temp = "t"+std::to_string(node->_type)+"_"+std::to_string(node->_index)+":"+std::to_string(node->_weight);
        std::string temp = "t"+std::to_string(static_cast<long long>(node->_type))+"_"+std::to_string(static_cast<long long>(node->_locX))+"_"+std::to_string(static_cast<long long>(node->_locY))+":"+std::to_string(static_cast<long long>(node->_weight));
        *lineage += temp;
        return lineage;
    }
    
    if (node->_children.size()==1) {
        lineage = getNWK(node->_children[0],lineage);
    }else{
        std::cout<<"(";
    *lineage += "(";
    for (int i=0; i<node->_children.size(); i++) {
        print(node->_children[i]);
        lineage = getNWK(node->_children[i],lineage);
        if (i!=node->_children.size()-1) {
            std::cout<<",";
            *lineage += ",";
        }
    }
    std::cout<<")"<<node->_index<<":"<<node->_weight;
        *lineage += ")";
    //*lineage += ")"+std::to_string(static_cast<long long>(node->_index))+":"+std::to_string(static_cast<long long>(node->_weight));
    }
    return lineage;
}

