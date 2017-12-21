//
//  Lattices.h
//  CancerModel_cleanup
//
//  Created by Yuantong Ding on 2/26/14.
//  Copyright (c) 2014 Yuantong Ding. All rights reserved.
//

#ifndef __CancerModel_cleanup__Lattices__
#define __CancerModel_cleanup__Lattices__

#include <iostream>
#include <vector>
#include <map>
#include "GenealogyNode.h"
#include "randgen.hpp"
#include "Tree.h"
#include <math.h>
#include "Topology.h"
#include <algorithm>



using namespace std;

class Lattices{
    std::vector<std::vector<GenealogyNode*>> cells;
    int _xStep;
    int _yStep;
    GenealogyNode * root_normal;
    GenealogyNode * root_tumor;
    ginkgo::RandomNumberGenerator& rng_;
    int _aliveCellCount;
    int _normalCellCount;
    double _emptyCellFittness;
    
    double _r;
    int _centerX;
    int _centerY;
    int _initialCellNumber;
    bool _hitBoundary;
    int _migrationEvent;
    int _tumorXmin;
    int _tumorYmin;
    int _tumorXmax;
    int _tumorYmax;
    Tree * _treeTumor;
    Tree * _treeNormal;
    std::vector<std::vector<CellIndexType>> _base;
    
public:
    Lattices(int xStep,int yStep, int initialCellNumber);
    
    Lattices(int xStep, int yStep);
    
    ~Lattices();
    
    void setRandomSeed(unsigned long seed);
    
    void printType();
    
    void printCell();
    
    void setEmptyCellFit(double fit);
        
    char decideFate(int x,int y,double t);
    
    bool migrate(int x,int y, double t);
    
    bool mutate(int x,int y,double t);
    
    bool replace(double fit, int n, int m);
    
    vector<int> emptyNeighbour(int x,int y,double t);
    
    vector<int> randomNeighbour(int x,int y,double t);
    
    vector<int> randomNeighbour2(int x,int y,double t);
    
    void proliferate(int x,int y, double t);
    
    vector<int> randomAliveCell();
    
    vector<int> randomSampleAreaLimit();
    
    vector<CellIndexType> getLineage(int x,int y);
    
    GenealogyNode * getCell(int x,int y);
        
    int normalCellNumber(){return _normalCellCount;};
    
    int aliveCellNumber(){return _aliveCellCount;};
    
    bool allNormalCell();
    
    void updateOxygen();
    
    void clear();
    
    int getCenterX(){return _centerX;};
    
    int getCenterY(){return _centerY;};
    
    double getSampleR(){return _r;};
    
    bool hitBoundary(){return _hitBoundary;};
    
    void sampling(int sampleSize, string filename,double time);
    
    vector<int> sampleTumor(int sampleSize, string samplefilename,double time, string seqFilename, string lineageFilename);
    
    void sampleGroup(int sampleSize, string samplefilename,double time, string seqFilename, string lineageFilename,string conFilename);
    
    void sampleGroup(vector<int>loc, int sampleSize, string samplefilename,double time, string seqFilename, string lineageFilename,string conFilename);
    
    void sampleSection(int sampleSize, string samplefilename,double time, string seqFilename, string lineageFilename);
    
    std::vector<std::vector<int>> sampleLayer(int sampleSize, string samplefilename,double time);
    
 //   void sampleLayer(int sampleSize, string samplefilename,double time, string seqFilename, string lineageFilename);
    
    vector<int> sampleTransect(int sampleSize, string samplefilename,double time, string seqFilename, string lineageFilename);
    
    void snapshot(string filename,double t);
        
    void getBase(string filename);
    
    std::vector<int> countCellNumber();
    
    void saveSampleLineage(std::vector<int> locX,std::vector<int> locY, string filename);
    
    void NJtree(vector<std::string> seqFileName);
    
//    void simulateSeq(std::vector<double> rate,Topology * t,int size);
    
private:
    Topology * nwk(int sampleSize,vector<int> X,vector<int> Y);
    Topology * nwk_noSingle(int sampleSize,vector<int> X,vector<int> Y);
    Topology *  nwk(int sampleSize,vector<int> X,vector<int> Y,std::vector<double> rate,string seqfilename,string lineageFilename);
    void simulateSeq(Topology * t, std::vector<double> rate, int length);
    void saveSeq(vector<int> X,vector<int> Y,Topology * t, string filename);
    
    void saveLineage(Topology *t, string filename);
    void saveLineage_noSingle(Topology *t, string filename);
    double adjustProliferationRate(int x, int y);
    std::string consensusSeq(vector<int> X,vector<int> Y, Topology * t);
    std::string consensusBase(std::vector<char>);
    void generateNewType(int current, int next);
    vector<int> sampleLoc(int min, int max, int sampleSize);

//    string randomBase(int n);
//    string pointMutate(string s);

};


#endif /* defined(__CancerModel_cleanup__Lattices__) */
