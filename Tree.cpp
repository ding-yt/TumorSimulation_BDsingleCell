//
//  Tree.cpp
//  CM
//
//  Created by Yuantong Ding on 12/14/14.
//  Copyright (c) 2014 Yuantong Ding. All rights reserved.
//

#include "Tree.h"


void Tree::setRoot(int name, int type){
    _root->_name = name;
    _root->_type = type;
    _root->_weight = 1;
    _root->_parent = NULL;
}

void Tree::addLineage(vector<CellIndexType> lineage){
    int start = lineage.size()-4;
    if (lineage[start]==_root->_name) {
        Node *current = _root;
        bool skip = false;
        for (int i=start-2; i>=0; i-=2) {
            for (int j=0; j<current->_children.size(); j++) {
                if (current->_children[j]->_name == lineage[i]) {
                    skip = true;
                    current = current->_children[j];
                    break;
                }
            }
            if (!skip) {
                Node *add;
                add->_name = lineage[i];
                add->_type = lineage[i+1];
                add->_weight = 1;
                add->_parent = current;
                current->_children.push_back(add);
            }
        }
    }
}

void Tree::compress(Node * node){
    Node *current = node;
    if (current->_children.size()==0) {
        return;
    }
    if (current->_children.size()==1 && current->_type == current->_children[0]->_type) {
        Node * p = current->_parent;
        Node * c = current->_children[0];
        current->_children[0]->_weight += current->_weight;
        current->_children[0]->_parent = current->_parent;
        for (int i=0; i<p->_children.size(); i++) {
            if (p->_children[i]->_name == current->_name) {
                p->_children[i] = current->_children[0];
                break;
            }
        }
        current->_parent = NULL;
        current->_children[0] = NULL;
        delete current;
        current = p;
    }else{
        for (int i=0; i<current->_children.size(); i++) {
            compress(current->_children[i]);
        }
    }
}

void Tree::compression(){
    compress(_root);
}


void Tree::print(Node * node){
    if (node->_children.size()==0) {
        cout<<node->_name<<":"<<node->_type<<":"<<node->_weight;
        return;
    }
    
    for (int i=0; i<node->_children.size(); i++){
    	if (i==0){
    		cout<<"(";
    	}
    	print(node->_children[i]);
    	cout <<", ";
    	if (i==node->_children.size()-1){
    		cout<<")";
    	}
    	
    }
    cout<<node->_name<<":"<<node->_type<<":"<<node->_weight;
}


void Tree::printTree(){
    print(_root);
}