//
//  SeqGen.h
//  CM
//
//  Created by Yuantong Ding on 2/3/15.
//  Copyright (c) 2015 Yuantong Ding. All rights reserved.
//

#ifndef __CM__SeqGen__
#define __CM__SeqGen__

#include <iostream>
#include <vector>

class SeqGen{
    std::vector<std::string> _seq;
    
public:
    SeqGen();
    void generateSeq(std::vector<double> rate);
    
};

#endif /* defined(__CM__SeqGen__) */
