/* theorem.h
   Copyright (C) 2009, Carlos Bedregal, all rights reserved.

   Implementation of Compressed Representation of Permutations: Runs & SRuns.

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

#ifndef THEOREM_H_INCLUDED
#define THEOREM_H_INCLUDED

#include "wavelettree.h"
#include "waveletnode.h"

using namespace std;

/** Base [abstract] class for the Hu-Tucker shaped Wavelet Tree [1] implementation.
 *
 *  [1] J. Barbay and G. Navarro, Compressed Representation of Permutations,
 *  and Applications.
 *
 *  @author Carlos Bedregal
 */

class Theorem{
    public:

    virtual ~Theorem(){};
    virtual uint length(){return len;}
    virtual WaveletTree<int> * tree() = 0;

    virtual uint pi(int i) = 0;
    virtual uint piInv(int i) = 0;

    /* saves the structure into files with prefix "fname" */
    virtual int save(char* fname) = 0;
    /* loads in memory a previously saved structure from file "fname" */
    virtual int load(char* fname) = 0;

    /* returns the number of bytes in memory */
    virtual int size() =0;
    /* returns the number of bits required by the bitsequences*/
    virtual unsigned int bitsRequired() =0;

    protected:
	/* size of permutation (aka size of tree's root) */
    uint len;
};

#endif // THEOREM_H_INCLUDED
