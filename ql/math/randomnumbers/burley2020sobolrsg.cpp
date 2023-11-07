/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2023 Peter Caspers

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include <ql/math/randomnumbers/burley2020sobolrsg.hpp>
#include <ql/math/randomnumbers/mt19937uniformrng.hpp>
#include <boost/functional/hash.hpp>

namespace QuantLib {

    Burley2020SobolRsg::Burley2020SobolRsg(Size dimensionality,
                                           unsigned long seed,
                                           SobolRsg::DirectionIntegers directionIntegers,
                                           unsigned long scrambleSeed)
    : dimensionality_(dimensionality), seed_(seed), directionIntegers_(directionIntegers),
      integerSequence_(dimensionality), sequence_(std::vector<Real>(dimensionality), 1.0) {
        reset();
        group4Seeds_.resize((dimensionality_ - 1) / 4 + 1);
        MersenneTwisterUniformRng mt(scrambleSeed);
        for (auto& s : group4Seeds_)
            s = static_cast<std::uint32_t>(mt.nextInt32());
    }

    void Burley2020SobolRsg::reset() const {
        sobolRsg_ = ext::make_shared<SobolRsg>(dimensionality_, seed_, directionIntegers_, false);
        nextSequenceCounter_ = 0;
    }

    const std::vector<std::uint_least32_t>& Burley2020SobolRsg::skipTo(std::uint32_t n) const {
        reset();
        for (Size k = 0; k < n + 1; ++k) {
            nextInt32Sequence();
        }
        return integerSequence_;
    }

    namespace {

        // for reverseBits() see http://graphics.stanford.edu/~seander/bithacks.html#BitReverseTable

        static const unsigned char bitReverseTable[] = {
            0,  128, 64, 192, 32, 160, 96,  224, 16, 144, 80, 208, 48, 176, 112, 240,
            8,  136, 72, 200, 40, 168, 104, 232, 24, 152, 88, 216, 56, 184, 120, 248,
            4,  132, 68, 196, 36, 164, 100, 228, 20, 148, 84, 212, 52, 180, 116, 244,
            12, 140, 76, 204, 44, 172, 108, 236, 28, 156, 92, 220, 60, 188, 124, 252,
            2,  130, 66, 194, 34, 162, 98,  226, 18, 146, 82, 210, 50, 178, 114, 242,
            10, 138, 74, 202, 42, 170, 106, 234, 26, 154, 90, 218, 58, 186, 122, 250,
            6,  134, 70, 198, 38, 166, 102, 230, 22, 150, 86, 214, 54, 182, 118, 246,
            14, 142, 78, 206, 46, 174, 110, 238, 30, 158, 94, 222, 62, 190, 126, 254,
            1,  129, 65, 193, 33, 161, 97,  225, 17, 145, 81, 209, 49, 177, 113, 241,
            9,  137, 73, 201, 41, 169, 105, 233, 25, 153, 89, 217, 57, 185, 121, 249,
            5,  133, 69, 197, 37, 165, 101, 229, 21, 149, 85, 213, 53, 181, 117, 245,
            13, 141, 77, 205, 45, 173, 109, 237, 29, 157, 93, 221, 61, 189, 125, 253,
            3,  131, 67, 195, 35, 163, 99,  227, 19, 147, 83, 211, 51, 179, 115, 243,
            11, 139, 75, 203, 43, 171, 107, 235, 27, 155, 91, 219, 59, 187, 123, 251,
            7,  135, 71, 199, 39, 167, 103, 231, 23, 151, 87, 215, 55, 183, 119, 247,
            15, 143, 79, 207, 47, 175, 111, 239, 31, 159, 95, 223, 63, 191, 127, 255};

        std::uint32_t reverseBits(std::uint32_t x) {
            std::uint32_t y;
            unsigned char* p = (unsigned char*)(&x);
            unsigned char* q = (unsigned char*)(&y);
            q[3] = bitReverseTable[p[0]];
            q[2] = bitReverseTable[p[1]];
            q[1] = bitReverseTable[p[2]];
            q[0] = bitReverseTable[p[3]];
            return y;
        }

        std::uint32_t laine_karras_permutation(std::uint32_t x, std::uint32_t seed) {
            x += seed;
            x ^= x * 0x6c50b47cu;
            x ^= x * 0xb82f1e52u;
            x ^= x * 0xc7afe638u;
            x ^= x * 0x8d22f6e6u;
            return x;
        }

        std::uint32_t nested_uniform_scramble(std::uint32_t x, std::uint32_t seed) {
            x = reverseBits(x);
            x = laine_karras_permutation(x, seed);
            x = reverseBits(x);
            return x;
        }

    }

    const std::vector<std::uint32_t>& Burley2020SobolRsg::nextInt32Sequence() const {
        const auto& seq =
            sobolRsg_->skipTo(nested_uniform_scramble(nextSequenceCounter_, group4Seeds_[0]));
        std::copy(seq.begin(), seq.end(), integerSequence_.begin());
        bool updating = true;
        Size i = 0, group = 0;
        do {
            Size seed = group4Seeds_[group++];
            for (Size g = 0; g < 4; ++g) {
                boost::hash_combine(seed, g);
                integerSequence_[i] = nested_uniform_scramble(integerSequence_[i], seed);
                if (++i >= dimensionality_) {
                    updating = false;
                }
            }
        } while (updating);
        ++nextSequenceCounter_;
        return integerSequence_;
    }

    const SobolRsg::sample_type& Burley2020SobolRsg::nextSequence() const {
        const std::vector<std::uint32_t>& v = nextInt32Sequence();
        // normalize to get a double in (0,1)
        for (Size k = 0; k < dimensionality_; ++k)
            sequence_.value[k] = static_cast<double>(v[k]) / 4294967296.0;
        return sequence_;
    }

}
