///////////////////////////////////////////////////////////////////////////////
//
// GINKGO Phylogeographical Evolution Simulator.
//
// Copyright 2009 Jeet Sukumaran and Mark T. Holder.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include "randgen.hpp"

namespace ginkgo {

///////////////////////////////////////////////////////////////////////////////
// RandomNumberGenerator

//! returns current seed
unsigned long RandomNumberGenerator::get_seed() const {
    return this->seed_;
}

//! seeds using given seed
void RandomNumberGenerator::set_seed(unsigned long seed) {
    this->seed_ = seed;
    srand(this->seed_);
}

//! returns a uniform random real between 0 and 1
float RandomNumberGenerator::uniform_01() {
    return static_cast<float>(rand())/static_cast<float>(RAND_MAX);
}

//! returns a uniform random integer between >= a and <= b
long RandomNumberGenerator::uniform_int(int a, int b) {
    return (rand() % (b-a+1)) + a;
}

//! Gaussian distribution with mean=1 and sd=0
//! from Knuth, The Art of Computer Programming, Sec 3.4.1, Algorithm P
float RandomNumberGenerator::standard_normal() {

    // since this method generates two variates at a time,
    // we store the second one to be returned on the next call
    static float stored_variate = 0;
    static bool return_stored=false;

    if (return_stored) {
        return_stored = false;
        return stored_variate;
    }

    float u1;
    float u2;
    float v1;
    float v2;
    float s = 1;

    while (s >= 1.0) {
        u1 = this->uniform_01();
        u2 = this->uniform_01();
        v1 = 2.0 * u1 - 1.0;
        v2 = 2.0 * u2 - 1.0;
        s = pow(v1, 2) + pow(v2, 2);
    }

    float polar = sqrt( (-2 * log(s)) / s );
    float x1 = v1 * polar;
    stored_variate = v2 * polar;

    return x1;
}

//! Gaussian with given mean and sd
float RandomNumberGenerator::normal(float mean, float sd) {
    return this->standard_normal() * sd + mean;
}

//! Poisson r.v. with given rate
unsigned int RandomNumberGenerator::poisson(float rate) {
    const int MAX_EXPECTATION = 64;
    if (rate > MAX_EXPECTATION) {
        float r = rate/2.0;
        return this->poisson(r) + this->poisson(r);
    }
    float L = exp(-1.0 * rate);
    float p = 1.0;
    float k = 0.0;
    while (p >= L) {
        k += 1.0;
        p *= this->uniform_01();
    }
    return static_cast<unsigned int>(k - 1.0);
}

//! Returns index according to probability given by weights.
unsigned int RandomNumberGenerator::weighted_index_choice(const std::vector<float>& weights) {
    float sum = 0.0;
    for (std::vector<float>::const_iterator i = weights.begin();
            i != weights.end();
            ++i) {
        sum += *i;
    }
    float rnd = this->uniform_01() * sum;
    unsigned int len = weights.size();
    for (unsigned int i = 0; i < len; ++i) {
        rnd -= weights[i];
        if (rnd < 0) {
            return i;
        }
    }
    return (len-1);
}

///////////////////////////////////////////////////////////////////////////////
// The global instance of the RNG

RandomNumberGenerator RandomNumberGenerator::instance_;

///////////////////////////////////////////////////////////////////////////////
// RandomPointer

RandomPointer::RandomPointer(RandomNumberGenerator& rng)
    : rng_(rng) { }

std::ptrdiff_t RandomPointer::operator() (std::ptrdiff_t max) {
    return static_cast<std::ptrdiff_t>(this->rng_.uniform_int(0, max-1));
}

} // ginkgo namespace
