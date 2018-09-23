//
//  Topology.h
//  NodeTester
//
//  Created by Yuantong Ding on 1/27/15.
//  Copyright (c) 2015 Yuantong Ding. All rights reserved.
//

#ifndef __NodeTester__Topology__
#define __NodeTester__Topology__

#include <iostream>
#include <vector>
#include "randgen.hpp"
typedef int           CellIndexType;

struct Nodetree{
    Nodetree * _parent;
    int _index;
    int _weight;
    int _type;
    int _locX;
    int _locY;
    std::string _name;
    std::string _sequence;
    std::vector<Nodetree *> _children;
    
    
public:
    Nodetree();
    ~Nodetree();
    Nodetree(int index, int type, Nodetree * parent);
    void addName(int type, int x,int y);
};

class Topology{
    
private:
    Nodetree * _root;
    ginkgo::RandomNumberGenerator& _rng;
    
    
public:
    Topology();
    ~Topology();
    void addRoot(int index, int type);
    void destoryTree(Nodetree * node);
    Nodetree * getParent();
    void setParent(Nodetree * parent);
    void setLocation(Nodetree * node,int x, int y){node->_locX = x; node->_locY = y;};
    Nodetree * addNode(Nodetree * parentNode, int index, int type);
    void addLineage(std::vector<CellIndexType> lineage);
    Nodetree * searchNode(int index);
    void compress();
    void compress_noSingleton();
    void printNWK();
    Nodetree * getRoot(){return _root;};
    void addChild(Nodetree * child);
    void setSeq(std::string s){_root->_sequence = s;};
    void simulateSeq(std::vector<double> rate, int length);
    void setRandomSeed(unsigned long seed){_rng.set_seed(seed);};
    std::vector<std::string> printSeq();
    std::string getSequence(int index);
    std::string getNWK();
    std::string getNWK_noSingle();

    
    
    
private:
    Topology(const Nodetree&);
    const Topology& operator=(const Topology&);
    void decreaseRef();
    void increaseRef();
    Nodetree * searchNode(int index, Nodetree * startNode);
    Nodetree * searchChildren(int index, Nodetree * parentNode);
    void print(Nodetree * node);
    void compress(Nodetree * node);
    void compress_noSingleton(Nodetree * node);
    std::string randomBase(int n);
    char pointMutate(char s);
    void mutateSeq(Nodetree * node, std::vector<double> rates, int * count);
    void chopSeq(Nodetree * node, int length);
    std::vector<std::string> getSeq(Nodetree * node,std::vector<std::string> seq);
    std::string getSeq(Nodetree * node,int index);
    std::string * getNWK(Nodetree * node, std::string *lineage);
    std::string * getNWK_noSingle(Nodetree * node, std::string *lineage);
};

#endif /* defined(__NodeTester__Topology__) */
