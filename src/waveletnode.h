/* waveletnode.h
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

#ifndef WAVELETNODE_H_INCLUDED
#define WAVELETNODE_H_INCLUDED

#include<iostream>
#include<math.h>
#include<static_bitsequence.h>

using namespace std;

/** Auxiliar class to handle nodes of a wavelet tree like structure.
 *
 *  @author Carlos Bedregal
 */

class WTNode{
    public:
    static_bitsequence* bitseq; //estructure for rank & select
    WTNode* children[2]; //array of children: 0=left, 1=right

    public:
    WTNode();
    WTNode(int s);
    ~WTNode();
    void print();
    void createBitseq(uint* bitmap, uint size);
    static static_bitsequence* bitseqCreator(uint* bitmap, uint size);
    int save(FILE * fp);
    int load(FILE * fp);
    int size();
};

WTNode::WTNode(){
    children[0]=children[1]=0;
}

WTNode::WTNode(int s){
    //size=s;
    #ifdef VERBOSE
        cout<<"size: "<<s<<", bitmap["<<uint_len(s,1)<<"]\n";
    #endif //VERBOSE
    children[0]=children[1]=0;
    bitseq=0;
}

WTNode::~WTNode(){
    if(bitseq) delete bitseq;
}

void WTNode::print(){
    cout<<"["<<bitseq->length()<<"]: ";
    for(uint q=0; q<bitseq->length(); q++)
        cout<<bitseq->access(q)<<" ";
    cout<<endl;
}

void WTNode::createBitseq(uint* bitmap, uint size){
    //bitseq = WTNode::bitseqCreator(bitmap,size);
    bitseq = new static_bitsequence_brw32(bitmap,size,FACTOR);
}

static_bitsequence* WTNode::bitseqCreator(uint* bitmap, uint size){
    switch(bitseqFlag){//bitseqFlag defined at runtime
        case RRR:
            return (new static_bitsequence_rrr02(bitmap,size));
        case RRRL:
            return (new static_bitsequence_rrr02_light(bitmap,size));
        default:
            return (new static_bitsequence_brw32(bitmap,size,FACTOR));
    }
}

int WTNode::save(FILE * fp){
    #ifdef DEBUG3
        cout<<this<<": bitseq: len "<<bitseq->length()<<", bytes "<<bitseq->size()<<endl;
    #endif //DEBUG3
    return bitseq->save(fp);
}

int WTNode::load(FILE * fp){
	//bitseq = static_bitsequence::load(fp);

	switch(bitseqFlag){//bitseqFlag defined at runtime
		case RRR:
			bitseq = static_bitsequence_rrr02::load(fp); break;
		case RRRL:
			bitseq = static_bitsequence_rrr02_light::load(fp); break;
		default:
			bitseq = static_bitsequence_brw32::load(fp); break;
	}

    if(bitseq){
        #ifdef DEBUG2
            cout<<this<<": bitseq: len "<<bitseq->length()<<", bytes "<<bitseq->size()<<endl;
        #endif //DEBUG2
        return 0;
    }
	#ifdef DEBUG2
		cout<<"no bitseq created!!!\n";
	#endif //DEBUG2

    return -1;
}

int WTNode::size(){
    return bitseq->size() + sizeof(WTNode) - sizeof(uint);
}

#endif // WAVELETNODE_H_INCLUDED
