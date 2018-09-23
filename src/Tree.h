//
//  Tree.h
//  CM
//
//  Created by Yuantong Ding on 12/14/14.
//  Copyright (c) 2014 Yuantong Ding. All rights reserved.
//

#ifndef __CM__Tree__
#define __CM__Tree__

#include <iostream>
#include <vector>
#include "GenealogyNode.h"

using namespace std;

struct Node{
    CellIndexType _name;
    int _type;
    int _weight;
    Node * _parent;
    vector<Node*> _children;
};

class Tree{
    
    
public:
    Node *_root;
    
    void setRoot(int name, int type);
    void addLineage(vector<CellIndexType> lineage);
    void compression();
    void printTree();
    
private:
    void compress(Node * node);
    void print(Node * node);
};

//class Treenode{
//    int name;
//    int weitht;
//    Node * parent;
//    vector<Node*> children;
//public:
//    Treenode();
//
//};


#endif /* defined(__CM__Tree__) */
