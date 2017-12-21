//
//  Lattices.cpp
//  CancerModel_cleanup
//
//  Created by Yuantong Ding on 2/26/14.
//  Copyright (c) 2014 Yuantong Ding. All rights reserved.
//

#include "Lattices.h"
#include <fstream>
#include <sstream>
#include "randgen.hpp"
#include "Topology.h"
#include <math.h>

#define PI 3.141592653

using namespace std;


Lattices::Lattices(int xStep, int yStep, int initialCellNumber):
root_normal(new GenealogyNode()),
root_tumor(new GenealogyNode()),
_xStep(xStep),
_yStep(yStep),
_centerX(xStep/2),
_centerY(yStep/2),
_r(sqrt(initialCellNumber/PI)),
_initialCellNumber(initialCellNumber),
_emptyCellFittness(0),
_hitBoundary(false),
_migrationEvent(0),
_tumorXmin(xStep/2),
_tumorYmin(yStep/2),
_tumorXmax(xStep/2),
_tumorYmax(yStep/2),
rng_(ginkgo::RandomNumberGenerator::get_instance())
{
    double r2 = initialCellNumber/3.14;
    _aliveCellCount = 0;
    _normalCellCount = 0;
    for (int i=0; i<_xStep; i++) {
        vector<GenealogyNode*> g;
        cells.push_back(g);
        for (int j=0; j<_yStep; j++) {
            GenealogyNode * node = new GenealogyNode();
            cells[i].push_back(node);
            double ii = i-_xStep/2;
            double jj = j-_yStep/2;
            if (ii*ii+jj*jj<=r2) {
                cells[i][j]->set_parent(root_normal);
                _aliveCellCount ++;
                _normalCellCount ++;
            }else{
                cells[i][j]->die(0);
            }
            
        }
    }
    cells[_centerX][_centerY]->set_parent(root_tumor);
    cells[_centerX][_centerY]->set_type(1);
    _normalCellCount --;
    _treeTumor->setRoot(1, 1);
    _treeNormal->setRoot(0, 0);
    cout<<"root_normal: "<<root_normal->get_cell_index()<<"\n";
    cout<<"root_tumor: "<<root_tumor->get_cell_index()<<"\n";
}

void Lattices::setEmptyCellFit(double fit){
    _emptyCellFittness = fit;
}

Lattices::Lattices(int xStep, int yStep):
_xStep(xStep),
_yStep(yStep),
_centerX(xStep/2),
_centerY(yStep/2),
_emptyCellFittness(0),
_r(0),
_tumorXmin(xStep/2),
_tumorYmin(yStep/2),
_tumorXmax(xStep/2),
_tumorYmax(yStep/2),
_hitBoundary(false),
_migrationEvent(0),
_aliveCellCount(xStep*yStep),
_normalCellCount(xStep*yStep-1),
_initialCellNumber(xStep*yStep),
rng_(ginkgo::RandomNumberGenerator::get_instance())
{
    root_normal = new GenealogyNode();
    root_tumor = new GenealogyNode();
    for (int i=0; i<_xStep; i++) {
        vector<GenealogyNode*> g;
        cells.push_back(g);
        for (int j=0; j<_yStep; j++) {
            GenealogyNode * node = new GenealogyNode();
            cells[i].push_back(node);
            cells[i][j]->set_parent(root_normal);
        }
    }
    root_tumor->set_type(1);
    cells[_centerX][_centerY]->set_parent(root_tumor);
    cells[_centerX][_centerY]->set_type(1);
    cout<<"root_normal: "<<root_normal->get_cell_index()<<" type "<<root_normal->get_type()<<"\n";
    cout<<"root_tumor: "<<root_tumor->get_cell_index()<<" type "<<root_tumor->get_type()<<"\n";
    cout<<"tumor at ("<<_centerX<<","<<_centerY<<") "<<cells[_centerX][_centerY]->get_cell_index()<<",type "<<cells[_centerX][_centerY]->get_type()<<"\n\n";
}

Lattices::~Lattices(){
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            delete cells[i][j];
            cells[i][j]=NULL;
        }
    }
    root_normal = NULL;
    root_tumor = NULL;
}

void Lattices::clear(){
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            delete cells[i][j];
            cells[i][j]=NULL;
        }
    }
    GenealogyNode::set_counter(0);
    root_normal = new GenealogyNode();
    root_tumor = new GenealogyNode();
    root_tumor->set_type(1);
//    root_normal->set_cell_index(0);
//    root_tumor->set_cell_index(1);
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            cells[i][j]= new GenealogyNode();
            cells[i][j]->set_parent(root_normal);
        }
    }
    cells[_centerX][_centerY]->set_parent(root_tumor);
    cells[_centerX][_centerY]->set_type(1);
    _aliveCellCount = _xStep*_yStep;
    _normalCellCount = _aliveCellCount-1;
    _tumorXmin = _xStep/2;
    _tumorYmin = _yStep/2;
    _tumorXmax = _xStep/2;
    _tumorYmax = _yStep/2;
    cout<<"root_normal: "<<root_normal->get_cell_index()<<" type "<<root_normal->get_type()<<"\n";
    cout<<"root_tumor: "<<root_tumor->get_cell_index()<<" type "<<root_tumor->get_type()<<"\n";

}

void Lattices::setRandomSeed(unsigned long seed){
    rng_.set_seed(seed);
}

void Lattices::printType(){
    std::cout<<"cell type:\n";
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            cout<<cells[i][j]->get_type()<<"\t";
        }
        cout <<"\n";
    }
}

void Lattices::printCell(){
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            cout<<cells[i][j]->get_cell_index()<<":"<<cells[i][j]->get_type()<<"\t";
        }
        cout <<"\n";
    }
}

GenealogyNode * Lattices::getCell(int x, int y){
    return cells[x][y];
}


/*
 * chose a random alive cell location from the lattices, return vector<int> loc, where loc[0] and loc[1]
 * are x and y coordinates respectively
 */
vector<int> Lattices::randomAliveCell(){
    vector<int> location;
    int rand = rng_.uniform_int(0, _xStep-1);
    location.push_back(rand);
    location.push_back(rng_.uniform_int(0, _yStep-1));
    while (!cells[location[0]][location[1]]->isAlive()) {
        location[0] = rng_.uniform_int(0, _xStep-1);
        location[1] = rng_.uniform_int(0, _yStep-1);
    }
    return location;
}

//////////////////??????????????????????????????
vector<int> Lattices::randomSampleAreaLimit(){
    vector<int> location;
    int original_r = sqrt(floor(_initialCellNumber/PI));
    int lowerX = _centerX-original_r;
    int upperX = _centerX+original_r;
    int lowerY = _centerY-original_r;
    int upperY = _centerY+original_r;
    int newR = _r;
    int rand = rng_.uniform_int(lowerX,upperX-1);
    
    location.push_back(rand);
    location.push_back(rng_.uniform_int(lowerY,upperY-1));
    double ii = location[0]-_centerX;
    double jj = location[1]-_centerY;
    double temp_r = ii*ii+jj*jj;
    
    //cout<<"random sample : "<<location[0]<<" "<<location[1]<<"original_r "<<original_r<<" _r "<<_r<<"\n";
    if (newR<original_r) {
        while (!cells[location[0]][location[1]]->isAlive() || temp_r>original_r*original_r) {
            location[0] = rng_.uniform_int(lowerX,upperX-1);
            location[1] = rng_.uniform_int(lowerY,upperY-1);
            ii = location[0]-_centerX;
            jj = location[1]-_centerY;
            temp_r = ii*ii+jj*jj;
        }
    }else{
        
        lowerX = _centerX-newR;
        upperX = _centerX+newR;
        lowerY = _centerY-newR;
        upperY = _centerY+newR;
        while (!cells[location[0]][location[1]]->isAlive() || temp_r>newR*newR) {
            location[0] = rng_.uniform_int(lowerX,upperX-1);
            location[1] = rng_.uniform_int(lowerY,upperY-1);
            ii = location[0]-_centerX;
            jj = location[1]-_centerY;
            temp_r = ii*ii+jj*jj;
        }
    }
    
    return location;
}

/*
 * print out and return cell lineage at a certain location, if cell is dead, print/return -1
 * lineage is a list of cell index of its parents and type
 * end with position of this cell
 * lineage: index of last child, type of last child, index of parents, type of parents,..., x, y
 */
vector<CellIndexType> Lattices::getLineage(int x,int y){
    vector<CellIndexType> indexList;
    if (!cells[x][y]->isAlive()) {
        //std::cout<<"\t*\n";
        indexList.push_back(-1);
        return indexList;
    }else{
        GenealogyNode* parent = cells[x][y];
        while (parent->get_parent()) {
            //            std::cout<<"\t"<<parent->get_cell_index();
            indexList.push_back(parent->get_cell_index());
            indexList.push_back(parent->get_type());
            parent = parent->get_parent();
        }
        //        std::cout<<"\t"<<parent->get_cell_index()<<"\n";
        indexList.push_back(parent->get_cell_index());
        indexList.push_back(parent->get_type());
        indexList.push_back(x);
        indexList.push_back(y);
        return indexList;
    }
}


/*
 * return the location of a random neighbour (4 neighbors) of a cell at certain location
 */
vector<int> Lattices::randomNeighbour(int x, int y,double t){
    vector<int> neighbors = {x-1,y,x,y+1,x,y-1,x+1,y};
    vector<int> neighbor;
    for (int i=0; i<8; i+=2) {
        if (neighbors[i]>=0 && neighbors[i]<_xStep && neighbors[i+1]>=0 && neighbors[i+1]<_yStep) {
            neighbor.push_back(neighbors[i]);
            neighbor.push_back(neighbors[i+1]);
        }
    }
    vector<int> randomNeighbourLocation;
    if (neighbor.size()>0) {
        int rand = rng_.uniform_int(0,neighbor.size()/2-1);        
        randomNeighbourLocation.push_back(neighbor[rand*2]);
        randomNeighbourLocation.push_back(neighbor[rand*2+1]);
        //cout<<"Neighbour size: "<<neighbor.size()<<"\t random choice: "<<rand<<"\n";
    }
    

    return randomNeighbourLocation;
}


/*
 * return the location of a random neighbour(8 neighbors) of a cell at certain location
 */
vector<int> Lattices::randomNeighbour2(int x, int y,double t){
    vector<int> neighbors = {x-1,y,x,y+1,x,y-1,x+1,y,x+1,y+1,x+1,y-1,x-1,y+1,x-1,y-1};
    vector<int> neighbor;
    for (int i=0; i<8; i+=2) {
        if (neighbors[i]>=0 && neighbors[i]<_xStep && neighbors[i+1]>=0 && neighbors[i+1]<_yStep) {
            neighbor.push_back(neighbors[i]);
            neighbor.push_back(neighbors[i+1]);
        }
    }
    int rand = rng_.uniform_int(0,neighbor.size()/2-1);
    vector<int> randomNeighbourLocation;
    randomNeighbourLocation.push_back(neighbor[rand*2]);
    randomNeighbourLocation.push_back(neighbor[rand*2+1]);
    //    cout<<"Neighbour size: "<<neighbour_count<<"\t random choice: "<<rand<<"\n";
    return randomNeighbourLocation;
}


/*
 * decide a cell at certain location will die('D') or Migrate('M') or proliferate ('P')
 */
char Lattices::decideFate(int x,int y,double t){
    if (!cells[x][y]->isAlive()) {
        return 'D';
    }
    double rand = rng_.uniform_01();
    double d = cells[x][y]->get_death_rate();
    double m = cells[x][y]->get_death_rate()+cells[x][y]->get_migration_rate();
//    double sum_rate = cells[x][y]->get_death_rate()+cells[x][y]->get_migration_rate()+cells[x][y]->get_proliferation_time()/cells[x][y]->get_maxP();
//    double d = cells[x][y]->get_death_rate()/sum_rate;
//    double p = (d+cells[x][y]->get_proliferation_time()/cells[x][y]->get_maxP())/sum_rate;
    //cout<<"cell at ("<<x<<", "<<y<<"), type "<<cells[x][y]->get_type()<<", d:"<<cells[x][y]->get_death_rate()<<", m"<<cells[x][y]->get_migration_rate()<<", p:"<<cells[x][y]->get_proliferation_time()/cells[x][y]->get_maxP()<<"\n";
    //cout<<sum_rate<<"\t"<<d<<"\t"<<p<<"\n";


    
    if (rand<d) {
        //cout<<rand<<" D\n";
        return 'D';
    }else if(rand<m){
        //cout<<rand<<" P\n";
        return 'M';
    }else{
        //cout<<d<<"\t"<<p<<"\tdie\t"<<rand<<"\n";
        //cout<<rand<<" M\n";
        return 'P';
    }
    
}

/*
 * For a cell at location (x,y), decide whether it will migrate
 */
bool Lattices::migrate(int x, int y, double t){
    if (!cells[x][y]->isAlive()) {
        return false;
    }
    double rand = rng_.uniform_01();
    //    cout<<"Migration rate: "<<cells.get(x, y)->get_migration_rate()<<"\t random: "<<rand<<"\n";//
    if (rand<cells[x][y]->get_migration_rate()) {
        cells[x][y]->die(t);
        _aliveCellCount --;
        if (cells[x][y]->get_type()==0) {
            _normalCellCount --;
        }
        return true;
    }else{
        return false;
    }
}

/*
 * For a cell at location (x,y) go through proliferation, decide whether it will mutate
 */
bool Lattices::mutate(int x, int y, double t){
    if (!cells[x][y]->isAlive()) {
        return false;
    }
    double rand = rng_.uniform_01();
    //    cout<<"Mutation rate: "<<cells.get(x, y)->get_mutation_rate()<<"\t random: "<<rand<<"\n";//
    if (rand<cells[x][y]->get_mutation_rate()) {
        // create new type if not exists
        int next = cells[x][y]->get_type()+1;
        //std::cout<<"mutate from "<<cells[x][y]->get_type()<<" to "<<next<<", type count "<<GenealogyNode::cellTypeCount()<<"\n";
        if (!cells[x][y]->existType(next)) {
            generateNewType(cells[x][y]->get_type(), next);
            std::cout<<"\nadd new type "<<next<<"\n";
            //GenealogyNode::show_allCellType();
        }
        return true;
    }else{
        return false;
    }
}

bool Lattices::replace(double fit, int n, int m){
    
    if (cells[n][m]->isAlive()) {
        double f = fit/(fit+cells[n][m]->get_fittness());
        double rand = rng_.uniform_01();
        //cout<<"fitness "<<fit<<" replace "<<cells[n][m]->get_fittness()<<", ratio "<<f<<", rand "<<rand<<"\n";
        if (rand<f) {
            //cout<<"\t sucess\n";
            return true;
        }else{
            return false;
        }
    }else{
        return true;
    }
}


/*
 * For a cell at location (x,y), return its dead neighbour location as x-y key-value pair
 * if there's no dead neighbour, map size is 0
 */
vector<int> Lattices::emptyNeighbour(int x,int y,double t){
    vector<int> emptyNeighbours;
    vector<int> neighbors = {x-1,y,x,y+1,x,y-1,x+1,y,x+1,y+1,x+1,y-1,x-1,y+1,x-1,y-1};
    vector<int> neighbor;
    for (int i=0; i<4; i+=2) {
        if (neighbors[i]>0 && neighbors[i]<_xStep && neighbors[i+1]>0 && neighbors[i+1]<_yStep && !cells[neighbor[i]][neighbor[i+1]]->isAlive()) {
            neighbor.push_back(neighbors[i]);
            neighbor.push_back(neighbors[i+1]);
        }
    }
    return emptyNeighbours;
}

/*
 * For a cell at location (x,y), if go through proliferation, decide whether its 2 offsprings mutate
 * The first offspring will occupy its random dead neighbour's space if available, or compete with a random alive neighbour
 * The second offspring will replace the parent's space
 */
void Lattices::proliferate(int x,int y, double t){
    char stage = decideFate(x, y, t);
    //cout<<"cell at ("<<x<<","<<y<<"), type "<<cells[x][y]->get_type()<<" stage "<<stage<<"\n";
    if (stage == 'P') {        
        double p = cells[x][y]->get_proliferation_time()/cells[x][y]->get_maxP();//###################
        //double p = cells[x][y]->get_proliferation_time();
//        if (cells[x][y]->get_type()>0) {
//            p = adjustProliferationRate(x, y);
//        }
       
        double rand = rng_.uniform_01();
        //cout<<"\tp "<<p<<" rand:"<<rand<<"\n";
        if (rand<p) {            
            // offspring mutate?
            GenealogyNode * temp_offsprint;
            bool mutated = mutate(x,y,t);//***********
            if (mutate(x, y, t)) {
                temp_offsprint = new GenealogyNode(cells[x][y]->get_type()+1,t);
            }else{
                temp_offsprint = new GenealogyNode(cells[x][y]->get_type(),t);
            }
            
            vector<int> neighbor = randomNeighbour(x, y, t);
           
            bool compete = replace(temp_offsprint->get_fittness(), neighbor[0], neighbor[1]);
            
            
            if (neighbor.size()>0) {
         
            int n_x = neighbor[0];
            int n_y = neighbor[1];
            
             //cout<<"\t p success, offspring type "<<temp_offsprint->get_type()<<"\n";
                //cout<<"\t replace cell at ("<<n_x<<","<<n_y<<"), type "<<cells[n_x][n_y]->get_type()<<"\n";
                
                
            if (compete) {
                temp_offsprint->set_parent(cells[x][y]);
                if ( n_x==0 || n_x==_xStep-1 || n_y==0 || n_y==_yStep-1) {
                    if (temp_offsprint->get_type()>0) {
                        _hitBoundary = true;
                    }
                }
                if (cells[n_x][n_y]->isAlive()) {
                    if (cells[n_x][n_y]->get_type()==0 && temp_offsprint->get_type()>0) {
                        _normalCellCount --;
                    }
                    if (cells[n_x][n_y]->get_type()>0 && temp_offsprint->get_type()==0) {
                        _normalCellCount ++;
                    }
                    cells[n_x][n_y]->die(t);
                }else{
                    if (temp_offsprint->get_type()==0) {
                        _normalCellCount ++;
                    }
                    _aliveCellCount ++;
                }
                
                cells[n_x][n_y] = temp_offsprint;
                if (cells[n_x][n_y]->get_type()>0) {
                    if (n_x<_tumorXmin) {
                        _tumorXmin = n_x;
                    }else if (n_x>_tumorXmax){
                        _tumorXmax = n_x;
                    }
                    if (n_y<_tumorYmin) {
                        _tumorYmin = n_y;
                    }else if (n_y > _tumorYmax){
                        _tumorYmax = n_y;
                    }
                }
                
            }else{
                delete temp_offsprint;
            }
            }
            
        
      
        // self mutate?
        //        cout<<"self mutate?\n";//
            //if (compete) {
         
        GenealogyNode * tempself;
        if (mutate(x, y, t)) {
            if (cells[x][y]->get_type()==0) {
                _normalCellCount --;
            }
            tempself = new GenealogyNode(cells[x][y]->get_type()+1,t);
          
        }else{
            tempself = new GenealogyNode(cells[x][y]->get_type(),t);
        }
        tempself->set_parent(cells[x][y]);
        cells[x][y]= tempself;
            
        }
        //}
        
    }else{  //'D' or 'M'
        //cout<<"die\n";
        if(stage=='M'){
            _migrationEvent ++;
        }
        if (cells[x][y]->isAlive()) {
            if (cells[x][y]->get_type()==0) {
                _normalCellCount --;
            }
            cells[x][y]->set_death_time(t);
            //cells[x][y]->die(t);
            _aliveCellCount --;
        }
    }
}


bool Lattices::allNormalCell(){
    if (_aliveCellCount == _normalCellCount) {
        return true;
    }else{
        return false;
    }
}

void Lattices::sampling(int sampleSize, string filename,double time){
    vector<int> locationX;
    vector<int> locationY;
    cout<<"samplefile: "<<filename<<"\n";
    for (int i=0; i<sampleSize; i++) {
        int x = rng_.uniform_int(0, _xStep-1);
        int y = rng_.uniform_int(0, _yStep-1);
        if (!cells[x][y]->isAlive()) {
            i--;
        }else{
            locationX.push_back(x);
            locationY.push_back(y);
        }
    }

    saveSampleLineage(locationX, locationY, filename);
    
}

vector<int> Lattices::sampleTumor(int sampleSize, string samplefilename,double time, string seqFilename, string lineageFilename){
    vector<int> locationX;
    vector<int> locationY;
    
    cout<<"\n sampling random ......\n";
    for (int i=0; i<sampleSize; i++) {
        int x = rng_.uniform_int(_tumorXmin, _tumorXmax);
        int y = rng_.uniform_int(_tumorYmin, _tumorYmax);
        if (!cells[x][y]->isAlive()) {
            i--;
        }else{
            locationX.push_back(x);
            locationY.push_back(y);
        }
    }
    std::vector<double> rate;
    for (int i=0; i<GenealogyNode::cellTypeCount(); i++) {
        rate.push_back(GenealogyNode::getMutationRate(i));
    }

    saveSampleLineage(locationX, locationY, samplefilename);
    
    nwk(sampleSize, locationX, locationY,rate,seqFilename,lineageFilename);
    
    vector<int> loc;
    for (int i=0; i<sampleSize; i++) {
        loc.push_back(locationX[i]);
        loc.push_back(locationY[i]);
    }
    return loc;
    
}

void Lattices::snapshot(string filename,double t){
    cout<<"snapshot: "<<filename<<"\n\n";
    ofstream snapshot_file;
    snapshot_file.open(filename.c_str());
    snapshot_file<<"cell at time "<<t<<"\n";
    snapshot_file<<"migration event count "<<_migrationEvent<<"\n";
    vector<int> typeCount(GenealogyNode::cellTypeCount());
    
//    int normalCount =0;
//    int type1count = 0;
//    int type2count = 0;
//    int type3count = 0;
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
//            if (cells[i][j]->get_type()==0) {
//                normalCount++;
//            }else if(cells[i][j]->get_type()==1){
//                type1count++;
//            }else if(cells[i][j]->get_type()==2){
//                type2count++;
//            }else if(cells[i][j]->get_type()==3){
//                type3count++;
//            }
            typeCount[cells[i][j]->get_type()] ++;
        }
    }
//    snapshot_file<<"type_0_count\t"<<normalCount<<"\n";
//    snapshot_file<<"type_1_count\t"<<type1count<<"\n";
//    snapshot_file<<"type_2_count\t"<<type2count<<"\n";
//    snapshot_file<<"type_3_count\t"<<type3count<<"\n";
    for (int tp = 0; tp<GenealogyNode::cellTypeCount(); tp++) {
        std::ostringstream s;
        s << tp;
        std::string count = "type_" + s.str() + "_count";
        snapshot_file<<count<<typeCount[tp]<<"\n";
    }
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            snapshot_file<<i<<"\t"<<j<<"\t";
            if (cells[i][j]->isAlive()) {
                snapshot_file<<cells[i][j]->get_type()<<"\n";
            }else{
                //snapshot_file<<"-1\n";
                snapshot_file<<"-"<<cells[i][j]->get_type()<<"\n"; //Modify: added 9/8/15
            }
        }
    }
    snapshot_file.close();
    
//    string sfilename = filename;
//    sfilename.append(".symbol");
//    ofstream pic;
//    pic.open(sfilename.c_str());
//    for (int i=0; i<_xStep; i++) {
//        for (int j=0; j<_yStep; j++) {
//            
//            if (cells[i][j]->isAlive()) {
//                pic<<cells[i][j]->get_type();
//            }else{
//                pic<<".";
//            }
//            
//        }
//        pic<<"\n";
//    }
//    pic.close();
    
    
    string spfilename = filename.append(".symbol");
    cout<<"snapshot: preparing"<<spfilename<<"\n";

    ofstream sp;
    sp.open(spfilename.c_str());
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            if (i<601 && i >300 && j>300 && j<601 ) {
                switch (cells[i][j]->get_type()) {
                    case 1:
                        sp<<"A";
                        break;
                    case 2:
                        sp<<"B";
                        break;
                    case 3:
                        sp<<"C";
                        break;
                    default:
                        sp<<"N";
                        break;
                }
            
            }else{
                if (cells[i][j]->isAlive()) {
                    sp<<cells[i][j]->get_type()<<"\t";
                }else{
                    sp<<".\t";
                }
            }
        }
        sp<<"\n";
    }
    sp.close();
    
    string lfilename = filename.append(".allLineage.txt");
    cout<<"snapshot: preparing"<<lfilename<<"\n\n";
    ofstream lineage;
    lineage.open(lfilename.c_str());
    lineage<<"snapshot at time "<<t<<"\n";
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            if (cells[i][j]->get_type()==1) {
                vector<CellIndexType> lin = getLineage(i, j);
                lineage<<"loc\t"<<i<<"\t"<<j<<"\t"<<cells[i][j]->get_type()<<"\n";
                for (int i=0; i<lin.size(); i++) {
                    lineage<<lin[i]<<"\t";
                }
                lineage<<"\n";
            }

        }
        lineage<<"\n";
    }
    lineage.close();
    cout<<"snapshot: "<<lfilename<<"\n\n";
}

void Lattices::sampleSection(int sampleSize, string samplefilename,double time, string seqFilename, string lineageFilename){
    vector<int> locationX;
    vector<int> locationY;
    cout<<"\n sampling section ...... \n";
    if (_tumorXmax-_tumorXmin>_tumorYmax-_tumorYmin) {
        int section = (_tumorXmax-_tumorXmin)/5;
        for (int i=0; i<5; i++) {
            for (int j=0; j<sampleSize/5; j++) {
                int x = rng_.uniform_int(_tumorXmin+i*section, _tumorXmin+(i+1)*section);
                int y = rng_.uniform_int(_tumorYmin, _tumorYmax);
                if (!cells[x][y]->isAlive()) {
                    i--;
                }else{
                    locationX.push_back(x);
                    locationY.push_back(y);
                }
            }
        }
    }else{
        int section = (_tumorYmax-_tumorYmin)/5;
        for (int i=0; i<5; i++) {
            for (int j=0; j<sampleSize/5; j++) {
                int x = rng_.uniform_int(_tumorXmin, _tumorXmax);
                int y = rng_.uniform_int(_tumorYmin+i*section, _tumorYmin+(i+1)*section);
                if (!cells[x][y]->isAlive()) {
                    i--;
                }else{
                    locationX.push_back(x);
                    locationY.push_back(y);
                }
            }
        }
    }
    
    saveSampleLineage(locationX, locationY, samplefilename);
    
    std::vector<double> rate;
    for (int i=0; i<GenealogyNode::cellTypeCount(); i++) {
        rate.push_back(GenealogyNode::getMutationRate(i));
    }
    nwk(sampleSize, locationX, locationY,rate,seqFilename,lineageFilename);
    

    
}

std::vector<std::vector<int>> Lattices::sampleLayer(int sampleSize, string samplefilename,double time){
    vector<int> locationX_core;
    vector<int> locationY_core;
    vector<int> locationX_boundary;
    vector<int> locationY_boundary;
    std::cout<<"\n sampling layer ...... \n";
    
    int r;
    int x_trans = _tumorXmin+(_tumorXmax-_tumorXmin)/2;
    int y_trans = _tumorYmin+(_tumorYmax-_tumorYmin)/2;
    if (_tumorXmax-_tumorXmin>_tumorYmax-_tumorYmin) {
        r = (_tumorYmax-_tumorYmin)/2;
    }else{
        r = (_tumorXmax-_tumorXmin)/2;
    }
    int r2 = r*r/2;
    r = (int)r/sqrt(r);
    
    for (int i=0; i<sampleSize; i++) {
        int x = rng_.uniform_int(x_trans-r,x_trans+r);
        int y = rng_.uniform_int(y_trans-r,y_trans+r);

        if (!cells[x][y]->isAlive()) {
            i--;
        }else{
            int tx = x-x_trans;
            int ty = y-y_trans;
            if (tx*tx+ty*ty<=r2) {
                locationX_core.push_back(x);
                locationY_core.push_back(y);
            }
        }
    }
    
    for (int i=0; i<sampleSize; i++) {
        int x = rng_.uniform_int(_tumorXmin, _tumorXmax);
        int y = rng_.uniform_int(_tumorYmin, _tumorYmax);
        if (!cells[x][y]->isAlive()) {
            i--;
        }else{
            int tx = x-x_trans;
            int ty = y-y_trans;
            if (tx*tx+ty*ty>r2) {
                locationX_boundary.push_back(x);
                locationY_boundary.push_back(y);
            }else{
                i--;
            }
        }
    }
   
    std::vector<double> rate;
    for (int i=0; i<GenealogyNode::cellTypeCount(); i++) {
        rate.push_back(GenealogyNode::getMutationRate(i));
    }
    
    string core_fileName = samplefilename +"_coreEnd.txt";
    string core_seqFilename = samplefilename+"_core"+"_seq"+"END.txt";
    string core_lineageFilename = samplefilename+"_core"+"_lineage"+"END.txt";

    saveSampleLineage(locationX_core, locationY_core, core_fileName);
    nwk(sampleSize, locationX_core, locationY_core,rate,core_seqFilename,core_lineageFilename);

    string boundary_fileName = samplefilename +"_boundaryEnd.txt";
    string boundary_seqFilename = samplefilename+"_boundary"+"_seq"+"END.txt";
    string boundary_lineageFilename = samplefilename+"_boundary"+"_lineage"+"END.txt";

    saveSampleLineage(locationX_boundary, locationY_boundary,  boundary_fileName);
    nwk(sampleSize, locationX_boundary, locationY_boundary,rate,boundary_seqFilename,boundary_lineageFilename);
    
    vector<int> loc_core;
    for (int i=0; i<sampleSize; i++) {
        loc_core.push_back(locationX_core[i]);
        loc_core.push_back(locationY_core[i]);
    }
    
    vector<int> loc_boundary;
    for (int i=0; i<sampleSize; i++) {
        loc_boundary.push_back(locationX_boundary[i]);
        loc_boundary.push_back(locationY_boundary[i]);
    }
    
    vector<vector<int>> loc;
    loc.push_back(loc_core);
    loc.push_back(loc_boundary);
    return loc;
    
}




void Lattices::sampleGroup(int sampleSize, string samplefilename,double time, string seqFilename, string lineageFilename, string confilename){
    vector<int> locationX;
    vector<int> locationY;
    vector<int> centerX;
    vector<int> centerY;
    cout<<"\n sampling group ...... \n";
    for (int i=0; i<sampleSize; i++) {
        int x = rng_.uniform_int(_tumorXmin, _tumorXmax);
        int y = rng_.uniform_int(_tumorYmin, _tumorYmax);
        if (!cells[x][y]->isAlive()) {
            i--;
        }else{
            centerX.push_back(x);
            centerY.push_back(y);
            
            for (int ii=x-2; ii<=x+2; ii++) {
                for (int jj=y-2; jj<=y+2;jj++) {
                    if (ii>0 && ii<_xStep && jj>0 &&jj<_yStep && cells[ii][jj]->isAlive()) {
                        locationX.push_back(ii);
                        locationY.push_back(jj);
                    }
                   
                }
            }

        }
    }
    std::vector<double> rate;
    for (int i=0; i<GenealogyNode::cellTypeCount(); i++) {
        rate.push_back(GenealogyNode::getMutationRate(i));
    }
    
    saveSampleLineage(locationX, locationY, samplefilename);
    
    Topology * t = nwk(locationX.size(), locationX, locationY,rate,seqFilename,lineageFilename);
    
    std::vector<string> consensus;
    ofstream conSeq;
    conSeq.open(confilename.c_str());
    
    for (int k=0; k<centerX.size(); k++) {
        int x = centerX[k];
        int y = centerY[k];
        std::vector<int> m;
        std::vector<int> n;
        conSeq<<">t"<<getCell(x, y)->get_type()<<"_"<<x<<"_"<<y<<"\n";

        for (int i=x-2; i<=x+2; i++) {
            for (int j=y-2; j<=y+2; j++) {
                if (i>0 && i<_xStep && j>0 &&j<_yStep && cells[i][j]->isAlive()) {
                    m.push_back(i);
                    n.push_back(j);
                }
                
            }
        }
        std::string seq = consensusSeq( m, n, t);
        consensus.push_back(seq);
        conSeq<<seq<<"\n";
    }
    conSeq<<">root\n";
    conSeq<<t->getRoot()->_sequence<<"\n";
    conSeq.close();
    
    std::cout<<"consensus seq file: "<<confilename<<"\n";
    
}

void Lattices::sampleGroup(vector<int> loc, int sampleSize, string samplefilename,double time, string seqFilename, string lineageFilename, string confilename){
    vector<int> locationX;
    vector<int> locationY;
    vector<int> centerX;
    vector<int> centerY;
    vector<int> subSampleN_count;
    vector<std::string> allSeqFiles;
    int subsample_leftsize = 3; // x-subsample_leftsize to x+subsample_rightsize, (3+1+3)^2 = 49 cells per sample max
    int subsample_rightsize = 3;
    
//    for (int i=0; i<sampleSize; i++) {
//        centerX.push_back(loc[2*i]);
//        centerY.push_back(loc[2*i+1]);
//    }
    
    std::vector<double> rate;
    for (int i=0; i<GenealogyNode::cellTypeCount(); i++) {
        rate.push_back(GenealogyNode::getMutationRate(i));
    }
    
    cout<<"\n sampling group ...... \n";
    unsigned long sampleN = loc.size()/2;
    for (int i=0; i<sampleN; i++) {
        int x = loc[2*i];
        int y = loc[2*i+1];
        if (!cells[x][y]->isAlive()) {
            i--;
        }else{
            subSampleN_count.push_back(0);
            centerX.push_back(x);
            centerY.push_back(y);
            
            for (int ii=x-subsample_leftsize; ii<=x+subsample_rightsize; ii++) {
                for (int jj=y-subsample_leftsize; jj<=y+subsample_rightsize;jj++) {
                    if (ii>0 && ii<_xStep && jj>0 &&jj<_yStep && cells[ii][jj]->isAlive()) {
                        locationX.push_back(ii);
                        locationY.push_back(jj);
                        subSampleN_count[i] ++;
                    }
                    
                }
            }
            
//            std::string groupSeqFile = seqFilename;
//            std::string name = "group_seq";
//            std::size_t found = groupSeqFile.find(name);
//            if (found!=std::string::npos) {
//                std::string newname = "group_seq_"+std::to_string(i);
//                groupSeqFile.replace(found, name.length(), newname);
//                Topology * t = nwk(locationX.size(), locationX, locationY,rate,groupSeqFile,lineageFilename);
//                allSeqFiles.push_back(groupSeqFile);
//            }
            
//            locationX.clear();
//            locationY.clear();
            
        }
    }
    
    

    
    saveSampleLineage(locationX, locationY, samplefilename);
    
    Topology * t = nwk(locationX.size(), locationX, locationY,rate,seqFilename,lineageFilename);
    
    //Modify: added 9/8/15
    int subSampleIndex = 0;
    for (int subi=0; subi<sampleN; subi++) {
        int x = loc[2*subi];
        int y = loc[2*subi+1];
        
        if (!cells[x][y]->isAlive()) {
            subi--;
        }else{
            std::string groupSeqFile = seqFilename;
            std::string name = "Group_seq";
            std::size_t found = groupSeqFile.find(name);
            vector<int> subSample_X;
            vector<int> subSample_Y;
            if (found!=std::string::npos) {
                std::string newname = "Group_seq_"+std::to_string(static_cast<long long>(subi));
                groupSeqFile.replace(found, name.length(), newname);
                for (int j=subSampleIndex; j<subSampleIndex + subSampleN_count[subi]; j++) {
                    subSample_X.push_back(locationX[j]);
                    subSample_Y.push_back(locationY[j]);
                }
                saveSeq(subSample_X, subSample_Y, t, groupSeqFile);
                allSeqFiles.push_back(groupSeqFile);
            }

            
        }
        subSampleIndex += subSampleN_count[subi];
    }
    //Modify: added 9/8/15
    
    std::vector<string> consensus;
    ofstream conSeq;
    conSeq.open(confilename.c_str());
    
    for (int k=0; k<centerX.size(); k++) {
        int x = centerX[k];
        int y = centerY[k];
        std::vector<int> m;
        std::vector<int> n;
        conSeq<<">t"<<getCell(x, y)->get_type()<<"_"<<x<<"_"<<y<<"\n";
        
        for (int i=x-subsample_leftsize; i<=x+subsample_rightsize; i++) {
            for (int j=y-subsample_leftsize; j<=y+subsample_rightsize; j++) {
                if (i>0 && i<_xStep && j>0 &&j<_yStep && cells[i][j]->isAlive()) {
                    m.push_back(i);
                    n.push_back(j);
                }
                
            }
        }
        std::string seq = consensusSeq( m, n, t);
        consensus.push_back(seq);
        conSeq<<seq<<"\n";
    }
    conSeq<<">root\n";
    conSeq<<t->getRoot()->_sequence<<"\n";
    conSeq.close();
    
    std::cout<<"consensus seq file: "<<confilename<<"\n";
    
    allSeqFiles.push_back(confilename);
    NJtree(allSeqFiles);
    
}

vector<int> Lattices::sampleTransect(int sampleSize, string samplefilename,double time, string seqFilename, string lineageFilename){

     cout<<"\n sampling section...... \n";
    
    vector<int> locationX;
    vector<int> locationY;
    vector<int> loc;
 
    int centerX = (_tumorXmax+_tumorXmin)/2;
    int centerY = (_tumorYmax+_tumorYmin)/2;
    if (_tumorXmax-_tumorXmin>_tumorYmax-_tumorYmin) {
        int range = (_tumorXmax-_tumorXmin)/2;
        int intvel = (int)floor(range/sampleSize);
        int i=0;
        for (int k=0; k<sampleSize; k++) {
            int x = centerX+i*intvel;
            int y = centerY;
            if (x>0 && x<_xStep) {
            if (!cells[x][y]->isAlive()) {
                k--;
            }else{
                locationX.push_back(x);
                locationY.push_back(y);
                loc.push_back(x);
                loc.push_back(y);
            }
            i++;
            }else{
                cout<<"!! section sampling number is "<<locationX.size()<<"\n";
                break;
            }
        }
    }else{
        int range = (_tumorYmax-_tumorYmin)/2;
        int intvel = (int)floor(range/sampleSize);
        int i=0;
        for (int k=0; k<sampleSize; k++) {
            int x = centerX;
            int y = centerY+i*intvel;
            if ( y>0 &&y<_yStep) {
            
                if (!cells[x][y]->isAlive()) {
                    k--;
                }else{
                    locationX.push_back(x);
                	locationY.push_back(y);
                    loc.push_back(x);
                    loc.push_back(y);
                }
                i++;
            }else{
                cout<<"!! section sampling number is "<<locationX.size()<<"\n";
                break;
            }
        }
    }
    
    std::vector<double> rate;
    for (int i=0; i<GenealogyNode::cellTypeCount(); i++) {
        rate.push_back(GenealogyNode::getMutationRate(i));
    }
    
    saveSampleLineage(locationX, locationY, samplefilename);
    
    nwk(sampleSize, locationX, locationY,rate,seqFilename,lineageFilename);
    
    return loc;
}

vector<int> Lattices::sampleLoc(int min, int max, int sampleSize){
    for (int k=0; k<sampleSize; k++) {
        int x = rng_.uniform_int(_tumorXmin, _centerX);
    }
}

Topology * Lattices::nwk(int sampleSize,vector<int> X,vector<int> Y,std::vector<double> rate,string seqfilename, string lineagefilename){
    Topology * t = nwk(sampleSize, X, Y);
//    Topology * t_noSingle = nwk_noSingle(sampleSize, X, Y);
    
    simulateSeq(t, rate, 3000);
    
    saveSeq(X, Y, t, seqfilename);
    
    saveLineage(t, lineagefilename);
    
//    string lineagefilename_noSingle = lineagefilename+"_noS.txt";
//    saveLineage_noSingle(t_noSingle, lineagefilename_noSingle);
    
    return t;
}

Topology * Lattices::nwk(int sampleSize,vector<int> locationX,vector<int> locationY){

    Topology * tAll = new Topology();

    tAll->addRoot(-1, -1);
    
    for (int s=0; s<locationX.size(); s++) {
      
        int x = locationX[s];
        int y = locationY[s];
        vector<CellIndexType> lineage = getLineage(x, y);
        vector<CellIndexType> baselineage;
        lineage.pop_back();
        lineage.pop_back();
        
        int ancestorIndex = lineage[lineage.size()-1-1-2];
        int ax = (int) (ancestorIndex-1)%_xStep;
        int ay = (ancestorIndex-1)%_xStep;
        int baseLoc = ax*_xStep+ay+1;
        baselineage = _base[baseLoc];
        
        vector<CellIndexType> combinedLineage = lineage;
        combinedLineage.insert(combinedLineage.end(), baselineage.begin(),baselineage.end());
        combinedLineage.push_back(x);
        combinedLineage.push_back(y);
        tAll->addLineage(combinedLineage);

    }
    
    tAll->compress();
//    std::cout<<"whole tree:\n";
//    tAll->printNWK();
    
    return tAll;
}


Topology * Lattices::nwk_noSingle(int sampleSize,vector<int> locationX,vector<int> locationY){
    
    Topology * tNo_single = new Topology();

    tNo_single->addRoot(-1, -1);
    
    for (int s=0; s<locationX.size(); s++) {
        
        int x = locationX[s];
        int y = locationY[s];
        vector<CellIndexType> lineage = getLineage(x, y);
        vector<CellIndexType> baselineage;
        lineage.pop_back();
        lineage.pop_back();
        
        int ancestorIndex = lineage[lineage.size()-1-1-2];
        int ax = (int) (ancestorIndex-1)%_xStep;
        int ay = (ancestorIndex-1)%_xStep;
        int baseLoc = ax*_xStep+ay+1;
        baselineage = _base[baseLoc];
        
        vector<CellIndexType> combinedLineage = lineage;
        combinedLineage.insert(combinedLineage.end(), baselineage.begin(),baselineage.end());
        combinedLineage.push_back(x);
        combinedLineage.push_back(y);
        tNo_single->addLineage(combinedLineage);
        
    }
    
    
    //print out lineage without singleton 9/23/15
    tNo_single->compress_noSingleton();
//    std::cout<<"whole tree no single:\n";
//    tNo_single->printNWK();
    
    return tNo_single;
}


void Lattices::simulateSeq(Topology * t, std::vector<double> rate,  int length){
    t->simulateSeq(rate, length);
}

void Lattices::saveLineage(Topology *tree, string filename){
    std::cout<<"lineage file: "<<filename<<"\n";
    ofstream lineage;
    lineage.open(filename.c_str());
    std::string s = tree->getNWK();
    lineage << s<<";";
    lineage.close();
}

void Lattices::saveLineage_noSingle(Topology *tree, string filename){
    std::cout<<"lineage file: "<<filename<<"\n";
    ofstream lineage;
    lineage.open(filename.c_str());
    std::string s = tree->getNWK_noSingle();
    lineage << s<<";";
    lineage.close();
}

void Lattices::saveSeq(vector<int> X,vector<int> Y,Topology * t, string filename){
    std::cout<<"seq file: "<<filename<<"\n";
    ofstream simSeq;
    simSeq.open(filename.c_str());
    assert(X.size()==Y.size());
    
    
    for (int s=0; s<X.size(); s++) {
        int sampleCount = 0;
        int x = X[s];
        int y = Y[s];
        vector<CellIndexType> lineage = getLineage(x, y);
        simSeq<<">t"<<getCell(x, y)->get_type()<<"_"<<x<<"_"<<y<<"\n";
        simSeq<< t->getSequence(lineage[0])<<"\n";
    }
    simSeq<<">root\n";
    simSeq<<t->getRoot()->_sequence<<"\n";
    simSeq.close();
}


void Lattices::getBase(std::string filename){
    
    //std::vector<std::vector<CellIndexType>> base;
    std::cout<<"read base from: "<<filename<<" ......\n";
    std::ifstream file;
    std::string line;
    file.open(filename.c_str());
    if (file.is_open()) {
        
        while (std::getline(file, line)) {
            std::vector<CellIndexType> baseLineage;
            std::stringstream ss(line);
            CellIndexType n;
            if (ss>>n){
                baseLineage.push_back(-n);
            while (ss >> n){
                baseLineage.push_back(-n);
            }
                baseLineage.pop_back();
                baseLineage.pop_back();
            _base.push_back(baseLineage);
            }
        }
        
        file.close();
    }else{
        std::cout << "Unable to open file "<<filename<<"\n";
    }
    
    std::cout <<"lineage base read in\t"<<_base.size()<<"\n";

//    Topology * baseTree = new Topology();
//    baseTree->addRoot(-1, -1);
//    for (int i=0; i<50; i++) {
//        baseTree->addLineage(_base[i]);
//    }
    
//    std::vector<double> r = {0.03,0.03,0,0};
//    baseTree->simulateSeq(r, 500);
//    
//    for (int i=0; i<50; i++) {
//        std::cout<<"search for :"<<_base[i][0]<<"\n";
//        std::cout<<baseTree->getSequence(_base[i][0])<<"\n";
//    }
//    
//    ofstream baseSeq;
//    string baseSeqfile = "/Users/dyt/Dropbox/cancerEvolution/script/CM/baselineageseq.txt";
//    baseSeq.open(baseSeqfile.c_str());
//    cout<<"basefile: "<<baseSeqfile<<"\n";
//    
//    for (int i=0; i<50; i++) {
//        baseSeq<<">"<<i<<"\n";
//        baseSeq<<baseTree->getSequence(_base[i][0])<<"\n";
//    }
//    baseSeq<<">root\n";
//    baseSeq<<baseTree->getSequence(baseTree->getRoot()->_index);
//    baseSeq.close();
    
//    baseTree->compress();
//    baseTree->printNWK();
//    
//    std::string l = baseTree->getNWK();
//    std::cout<<"\n\n"<<l<<"\n\n";
    
}


std::vector<int> Lattices::countCellNumber(){
    int typeN = GenealogyNode::cellTypeCount();
    std::vector<int> counts;
    for (int i = 0; i<typeN; i++) {
        counts.push_back(0);
    }
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            int tempType = cells[i][j]->get_type();
            counts[tempType]++;

        }
    }
    return counts;
    
}

void Lattices::saveSampleLineage(std::vector<int> locX,std::vector<int> locY, string filename){
    
    if (locX.size() != locY.size()) {
        std::cout<<"warning!!!!! sample location unequal\n";
    }else{
    ofstream sample_file;

    sample_file.open(filename.c_str());
    cout<<"samplefile: "<<filename<<"\n";
    
    sample_file<<"sample at time "<<time<<"\n";
    
    std::vector<int> cellTypeCount = countCellNumber();
    for (int i=0; i<GenealogyNode::cellTypeCount(); i++) {
        std::string prefix = "type_";
        std::ostringstream s;
        s << i;
        std::string proliferation = prefix + s.str() + "_count";
        sample_file<<proliferation<<"\t"<<cellTypeCount[i]<<"\n";
    }
    
    for (int s=0; s<locX.size(); s++) {
        int sampleCount = 0;
        int x = locX[s];
        int y = locY[s];
        vector<CellIndexType> lineage = getLineage(x, y);
        sample_file<<"sample"<<s<<"\n";
        for (int i=0; i<lineage.size(); i++) {
             sample_file<<lineage[i]<<"\t";
        }
        sample_file<<"\n";
        
    }
     sample_file.close();
    }
    

}

double Lattices::adjustProliferationRate(int x, int y){
    int xMin = _xStep/2;
    int xMax = _xStep/2;
    int yMin = _yStep/2;
    int yMax = _yStep/2;
    double p = cells[x][y]->get_proliferation_time()/cells[x][y]->get_maxP();
    
    for (int i=0; i<_xStep; i++) {
        for (int j=0; j<_yStep; j++) {
            if (cells[i][j]->get_type()> 0) {
                if (i<xMin) {
                    xMin = i;
                }else if (i>xMax){
                    xMax = i;
                }
                if (j<yMin) {
                    yMin = j;
                }else if (j > yMax){
                    yMax = j;
                }
            }
        }
    }
    
    double tumor_centerX = (xMin + xMax)/2;
    double tumor_centerY = (yMin + yMax)/2;
    double tumor_rangeX = - xMin + xMax;
    double tumor_rangeY = - yMin + yMax;
    double tumor_r;
    if (tumor_rangeX > tumor_rangeY) {
        tumor_r = tumor_rangeY;
    }else{
        tumor_r = tumor_rangeX;
    }
    double s_x = (x-tumor_centerX);
    double s_y = (y-tumor_centerY);
    double sample_r =sqrt(s_x*s_x + s_y*s_y);
    
    //std::cout<<"tumor ("<<tumor_centerX<<","<<tumor_centerY<<") "<<tumor_r<<"\n";
    //std::cout<<"sample ("<<x<<","<<y<<") "<<sample_r<<"\tp"<<p<<"\t";
    
    if (tumor_r>_xStep/4) {
        double ratio = sample_r/tumor_r;
        if (ratio<1) {
            p = (sample_r/tumor_r) * p;
        }
    }

    return p;
    
}

std::string Lattices::consensusSeq(vector<int> X,vector<int> Y,Topology * t){

    assert(X.size()==Y.size());
    std::vector<string> seq;
    for (int s=0; s<X.size(); s++) {
        int x = X[s];
        int y = Y[s];
        seq.push_back(t->getSequence(cells[x][y]->get_cell_index())) ;
    }
    assert(seq.size()>0);
    unsigned long length = seq[0].length();
    std::string consensus = "";
    for (int i=0; i<length; i++) {
        std::vector<char> base;
        for (int j=0; j<seq.size(); j++) {
            base.push_back(seq[j].at(i));
        }
        consensus.append(consensusBase(base));
    }

    return consensus;
}

std::string Lattices::consensusBase(std::vector<char> base){
    int A_count = 0;
    int T_count = 0;
    int C_count = 0;
    int G_count = 0;
    int other = 0;
    int max = 0;
    std::string consensus;
    
    for (int i=0; i<base.size();i++) {
        if (base[i] == 'A') {
            A_count ++;
        }else if (base[i] == 'T'){
            T_count ++;
        }else if (base[i] == 'C'){
            C_count ++;
        }else if (base[i] == 'G'){
            G_count ++;
        }else{
            other++;
        }
    }
    
    if (A_count>max) {
        consensus = 'A';
        max = A_count;
    }
    if (T_count>max) {
        consensus = 'T';
        max = T_count;
    }
    if (C_count>max) {
        consensus = 'C';
        max = C_count;
    }
    if (G_count>max) {
        consensus = 'G';
        max = G_count;
    }
    
    return consensus;
    
}

void Lattices::generateNewType(int current, int next){
    int name = next;
    double proliferation_time = GenealogyNode::getProliferationRate(current);
    double death_rate = GenealogyNode::getDeathRate(current);
    double transition_rate = GenealogyNode::getTransitionRate(current);
    double migration_rate = GenealogyNode::getMigrationRate(current);
    double fittness = GenealogyNode::getFittness(current);
    double point_mutation_rate = GenealogyNode::getMutationRate(current);
    double mean;
    double sd;
    double r;
    
    double pratio = GenealogyNode::getProliferationRate(2)/GenealogyNode::getProliferationRate(1)-1;
    double fratio = (GenealogyNode::getFittness(2)-1)/(GenealogyNode::getFittness(1)-1)-1;
    
    long choice = rng_.uniform_int(1, 5);
        
  
//            mean = GenealogyNode::getProliferationRate(current);
//            sd = mean/10;
//            r = rng_.normal(mean, sd);
//            while (r<mean) {
//                r = rng_.normal(mean, sd);
//            }
            proliferation_time = GenealogyNode::getProliferationRate(current)*(1+pratio);

//            mean = GenealogyNode::getDeathRate(current);
//            sd = mean/10;
//            r = rng_.normal(mean, sd);
//            while (r<0) {
//                r = rng_.normal(mean, sd);
//            }
            death_rate = GenealogyNode::getDeathRate(current);
   
//            mean = GenealogyNode::getTransitionRate(current);
//            sd = mean/10;
//            r = rng_.normal(mean, sd);
//            while (r<0) {
//                r = rng_.normal(mean, sd);
//            }
            transition_rate =  GenealogyNode::getTransitionRate(current);;
  
//            mean = GenealogyNode::getMigrationRate(current);
//            sd = mean/10;
//            r = rng_.normal(mean, sd);
//            while (r<0) {
//                r = rng_.normal(mean, sd);
//            }
            migration_rate = GenealogyNode::getMigrationRate(current);
  
//            mean = GenealogyNode::getFittness(current);
//            sd = mean/10;
//            r = rng_.normal(mean, sd);
//            while (r<1 || r>mean) {
//                r = rng_.normal(mean, sd);
//            }
            fittness = 1+(GenealogyNode::getFittness(current)-1)*(1+fratio);

    
    GenealogyNode::addType(name, point_mutation_rate, proliferation_time, death_rate, migration_rate, fittness, transition_rate);
    
  
}

void Lattices::NJtree(vector<std::string> seqFileName){
    ofstream R_file;
    std::string filename = "/Users/dyt/Dropbox/cancerEvolution/script/CM/NJtree.R";
    R_file.open(filename.c_str());
    
    R_file<<"library(\"ape\", lib.loc=\"/Library/Frameworks/R.framework/Versions/3.1/Resources/library\")\n";
    for (int i=0; i<seqFileName.size(); i++) {
        std::string treeFile = seqFileName[i];
        treeFile.replace(treeFile.find("txt"), 3, "nwk");
        R_file<<"seq <- read.dna(\""<<seqFileName[i]<<"\",format=\"fasta\")\n";
        R_file<<"dist <- dist.dna(seq,model=\"JC\")\n";
        R_file<<"tree <- nj(dist)\n";
        R_file<<"tree <- root(tree,\"root\")\n";
        R_file<<"write.tree(tree,file=\""<<treeFile<<"\")\n\n";
    }
    R_file.close();
    std::cout<<"R file "<<filename<<" ready\n";
    
    
}


