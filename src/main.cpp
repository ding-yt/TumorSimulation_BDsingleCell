//
//  main.cpp
//  CM
//
//  Created by Yuantong Ding on 12/12/14.
//  Copyright (c) 2014 Yuantong Ding. All rights reserved.
//

#include <iostream>
#include "ParFile.h"
#include "Lattices.h"
#include "GenealogyNode.h"
#include "Topology.h"

int main(int argc, const char * argv[])
{
    
    std::string setting_file, output,base_lineage_file;
    std::map<std::string, double> parameters;
    int fail_count = 0;
    long random_seed;
    bool repeat_simulation = true;
    bool first_migration_event = true;
    map<double, vector<int>> migration_event;
    int migration_count = 0;
    int fileIndex = 0;
    ParFile p;
    
    clock_t start_time = clock();
    
    for (int i=1; i<argc; i++) {
        if (i+1 != argc) {
            if (strcmp(argv[i], "-i") == 0){
                setting_file = argv[i + 1];
            }else if (strcmp(argv[i], "-o") == 0){
                output = argv[i+1];
            }else if (strcmp(argv[i], "-s") == 0){
                random_seed = std::atoi(argv[i+1]);
            }else if (strcmp(argv[i], "-b") == 0){
                base_lineage_file = argv[i+1];
            }
        }
        
    }
    
    output.erase(std::remove(output.begin(), output.end(), '\n'), output.end());
    std::cout<<"output is "<<output<<"\n";
    
    parameters = p.get_parameters(setting_file);
    p.show();
    
    GenealogyNode::set_allCellType(parameters);
    GenealogyNode::show_allCellType();
    GenealogyNode::set_counter(0);
    
    double time_max = parameters["time_max"];
    int xStep = parameters["xStep"];
    int yStep = parameters["yStep"];
    double sample_time_interval = parameters["sample_time_interval"];
    int sample_size = (int)parameters["sample_size"];
    double time = 0;   
    
    Lattices space(xStep,yStep);
    space.setRandomSeed(random_seed);
    //space.getBase("/Users/dyt/Dropbox/cancerEvolution/script/CM/par_baseEND.txt.symbol.allLineage.txt");
    //space.getBase("/Users/dyt/Dropbox/cancerEvolution/script/CM/par_baseTest.txt");
    space.getBase(base_lineage_file);
   
    while (repeat_simulation) {
        space.proliferate(xStep/2, yStep/2, 1);
        for (time=2; time<time_max; time++) {
            //Boundary condition check
            if (space.hitBoundary()) {
                repeat_simulation = false;
                cout<<"hit boundary! time:"<<time<<"\n";
                goto output;
            }
            
            //all normal cell check
            if (space.allNormalCell()) {
                repeat_simulation = true;
                fail_count ++;
                cout<<"faill attemp: allNormalCell\t"<<fail_count<<"\ttime:"<<time<<"\taliveCell:"<<space.aliveCellNumber()<<"\n";
                //string snapshot_fileName = output+to_string(fail_count)+".txt";
                //space.snapshot(snapshot_fileName, time);
                break;
            }else{
                repeat_simulation = false;
            }
            
            //pick a random cell to proliferate
            vector<int> loc = space.randomAliveCell();
            space.proliferate(loc[0], loc[1], time);
            
            //sampling at certain time interval
            if ((int)time%(int)sample_time_interval == 0 && (!repeat_simulation)) {
            //if ((space.aliveCellNumber() - space.normalCellNumber())%2000 == 0 && (!repeat_simulation)) {
                // prepare output file
                std::string filename = output;
                std::ostringstream s;
                std::ostringstream intervel;
                intervel << sample_time_interval;
                std::string new_filename = filename + s.str()+".txt";
                fileIndex ++;
                s<<fileIndex;
                new_filename = filename + s.str()+"_i"+intervel.str()+".txt";
                               
                //sample file
                string sample_fileName = output+"_sample"+s.str()+".txt";
                string seq_fileName = output+"_seq"+s.str()+".txt";
                string lineage_fileName = output+"_lineage"+s.str()+".txt";

                space.snapshot(new_filename, time);
                space.sampleTumor(sample_size,sample_fileName,time,seq_fileName,lineage_fileName);
            }
            
        }
        space.clear();
        GenealogyNode::set_allCellType(parameters);
    }
    
output:
    // prepare output file
    // sampling after hit boundary
    std::string filename = output;
    cout<<"time"<<time<<"\n";
    
    string snapshot_fileName = output+"END.txt";
    
    string sample_fileName = output+"_rand"+"END.txt";
    string seqFilename = output+"_rand_seq"+"END.txt";
    string lineageFilename = output+"_rand_lineage"+"END.txt";
    
    string sampleTransect_fileName = output +"_secEnd.txt";
    string tran_seqFilename = output+"_sec"+"_seq"+"END.txt";
    string tran_lineageFilename = output+"_sec"+"_lineage"+"END.txt";
    
    string sampleGroup_fileName = output +"_randGroupEnd.txt";
    string group_seqFilename = output+"_randGroup"+"_seq"+"END.txt";
    string group_lineageFilename = output+"_randGroup"+"_lineage"+"END.txt";
    string group_consensus = output+"_randGroup"+"_conSeq"+"END.txt";
    
    string tranGroup_fileName = output +"_secGroupEnd.txt";
    string trangroup_seqFilename = output+"_secGroup"+"_seq"+"END.txt";
    string trangroup_lineageFilename = output+"_secGroup"+"_lineage"+"END.txt";
    string trangroup_consensus = output+"_secGroup"+"_conSeq"+"END.txt";
    
    string coreGroup_fileName = output +"_coreGroupEnd.txt";
    string coregroup_seqFilename = output+"_coreGroup"+"_seq"+"END.txt";
    string coregroup_lineageFilename = output+"_coreGroup"+"_lineage"+"END.txt";
    string coregroup_consensus = output+"_coreGroup"+"_conSeq"+"END.txt";
    
    string boundaryGroup_fileName = output +"_boundaryGroupEnd.txt";
    string boundarygroup_seqFilename = output+"_boundaryGroup"+"_seq"+"END.txt";
    string boundarygroup_lineageFilename = output+"_boundaryGroup"+"_lineage"+"END.txt";
    string boundarygroup_consensus = output+"_boundaryGroup"+"_conSeq"+"END.txt";
    
//    string sampleSection_fileName = output +"_sectionEnd.txt";
//    string sec_seqFilename = output+"_section"+"_seq"+"END.txt";
//    string sec_lineageFilename = output+"_section"+"_lineage"+"END.txt";
    
//    string sampleLayer_fileName = output +"_layerEnd.txt";
//    string Layer_seqFilename = output+"_layer"+"_seq"+"END.txt";
//    string Layer_lineageFilename = output+"_layer"+"_lineage"+"END.txt";



    GenealogyNode::show_allCellType();
    
    
    space.snapshot(snapshot_fileName, time);
    
    vector<int> location = space.sampleTumor(sample_size,sample_fileName,time,seqFilename,lineageFilename);
    space.sampleGroup(location,sample_size, sampleGroup_fileName, time, group_seqFilename, group_lineageFilename, group_consensus);
    
    vector<int> locationSection = space.sampleTransect(sample_size, sampleTransect_fileName, time, tran_seqFilename, tran_lineageFilename);
    space.sampleGroup(locationSection,sample_size, tranGroup_fileName, time, trangroup_seqFilename, trangroup_lineageFilename, trangroup_consensus);
    
    vector<vector<int>> locationLayer = space.sampleLayer(sample_size, output, time);
    space.sampleGroup(locationLayer[0],sample_size, coreGroup_fileName, time, coregroup_seqFilename, coregroup_lineageFilename, coregroup_consensus);
    space.sampleGroup(locationLayer[1],sample_size, boundaryGroup_fileName, time, boundarygroup_seqFilename, boundarygroup_lineageFilename, boundarygroup_consensus);
    
    
    clock_t temp_time = clock();
    double sec = (temp_time-start_time)/CLOCKS_PER_SEC;
    int min = (int)sec/60;
    sec = sec - min*60;
    cout << "\nRunning time: " << min <<"min"<<sec<<"sec\n";
    
    
    return 0;
    
}

