/* theorem2.h
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

#ifndef THEOREM2_H_INCLUDED
#define THEOREM2_H_INCLUDED

#include"theorem.h"
#include"theorem1.h"

/** Implementation of Compressed Data Structure for Permutations based on 
 *	Wavelet Tree (Hu-Tucker) using strict ascending sub-sequences (SRuns).
 *  (Theorem2 or TH2 for practical use)
 *
 *  [1] J. Barbay and G. Navarro, Compressed Representation of Permutations,
 *  and Applications.
 *
 *  @author Carlos Bedregal
 */

class Theorem2:public Theorem{
    public:
    Theorem1 *th1;
    static_bitsequence* bitseqR;
    static_bitsequence* bitseqRinv;

    public:
    Theorem2();
    Theorem2(Permutation<int> *p);
    virtual ~Theorem2();

    WaveletTree<int> * tree();

    uint pi(int i);
    uint piInv(int i);

    int save (char* fname);
    int load (char* fname);

    int size();
    unsigned int bitsRequired();
};

Theorem2::Theorem2(){
    th1=0;
    bitseqR=0;
    bitseqRinv=0;
}

Theorem2::Theorem2(Permutation<int> *p){
    assert(p!=0);
    assert(p->len>0);

    #ifdef PRINT
        cout<<"+ Building Th2: "<<p->len<<" elements\n";
    #endif //PRINT

    //Permutation<int> *p=new Permutation<int>(array,n);
    len=p->len;

    int i, pos=0, szBMap, len_;
    unsigned int ii;
    uint* R;
    uint* Rinv;

    //Set bitmap R
    #ifdef PRINT
        cout<<"\t- Th2: building bitmap R\n";
    #endif //PRINT
	p->findSRuns();
	szBMap = uint_len(len,1);
    R=new uint[szBMap];
    for(i=0; i<szBMap; R[i++]=0);
    bitset(R,pos);
    for(i=0; i<p->tau-1; i++){
        pos=pos+p->SRuns[i];
        bitset(R,pos);
    }
	delete[]p->SRuns;
	len_ = p->tau;
	p->tau = 0;

    //Create inverse permutation of array
    int *arrayInv=new int[len];
    for(ii=0;ii<len;ii++)
        arrayInv[p->array[ii]]=ii;

    //Set bitmap Rinv
    #ifdef PRINT
        cout<<"\t- Th2: building bitmap Rinv\n";
    #endif //PRINT
    Rinv=new uint[szBMap];
    for(i=0; i<szBMap; Rinv[i++]=0);
    for(ii=0; ii<len; ii++)
        if(bitget(R,arrayInv[ii]))
            bitset(Rinv,ii);

    delete[]arrayInv;

    bitseqR = WTNode::bitseqCreator(R,len);
    delete[]R;
    bitseqRinv = WTNode::bitseqCreator(Rinv,len);
    delete[]Rinv;

    //create permutation' of size [tau]
    int *array_ = new int[len_];
    for(i=0; i<len_; i++){
        array_[i]=bitseqRinv->rank1(p->array[bitseqR->select1(i+1)])-1;
    }

	p=0;

    Permutation<int> *p_=new Permutation<int>(array_,len_);
	p_->findRuns();
    th1 = new Theorem1(p_);

    delete[]array_;
    delete p_;
}

Theorem2::~Theorem2(){
    delete bitseqR;
    delete bitseqRinv;
    if(th1!=0) delete th1;
}

WaveletTree<int>* Theorem2::tree(){
    return th1->wt;
}

uint Theorem2::pi(int i){
    int i_, j_;
    i_ = bitseqR->rank1(i)-1;
    j_ = th1->pi(i_);

    i_ = bitseqR->select1(i_+1);
    j_ = bitseqRinv->select1(j_+1);

    return j_ + i - i_;
}

uint Theorem2::piInv(int i){
    int i_, j_;
    i_ = bitseqRinv->rank1(i)-1;
    j_ = th1->piInv(i_);

    i_ = bitseqRinv->select1(i_+1);
    j_ = bitseqR->select1(j_+1);

    return j_ + i - i_;
}

/* saves structure TH2's bitmaps into files with prefix "fname"
 * - first: th1 structure
 * - second: bitsequence R
 * - third: bitsequence Rinv
 */
int Theorem2::save (char* fname){
	int ret = th1->save(fname);
	FILE * output;
    output = fopen(fname,"ab");
    if(bitseqR->save(output)!=0) return -1;
    if(bitseqRinv->save(output)!=0) return -1;
    fclose(output);
    return ret;
}

/* loads in memory a previously saved TH2 structure from file "fname"
 * - fname contains both the th1 structure and the R and Rinv bitmaps
 * - fname.idx: contains the three shape
 */
int Theorem2::load (char* fname){
	char fname2[128];
	strcpy(fname2,fname);
	strcat(fname2,".idx");

    FILE * input;
    input = fopen(fname,"rb");
	FILE * hierarchy;
	hierarchy = fopen(fname2,"rb");

    th1 = new Theorem1();
    int ret = th1->loadWT(input,hierarchy);
	if(ret==-1){
		cout<<"@Theorem1::load()\n";
		return -1;
	}

	switch(bitseqFlag){//bitseqFlag defined at runtime
		case RRR:
			bitseqR = static_bitsequence_rrr02::load(input);
			bitseqRinv = static_bitsequence_rrr02::load(input);
			break;
		case RRRL:
			bitseqR = static_bitsequence_rrr02_light::load(input);
			bitseqRinv = static_bitsequence_rrr02_light::load(input);
			break;
		default:
			bitseqR = static_bitsequence_brw32::load(input);
			bitseqRinv = static_bitsequence_brw32::load(input);
			break;
	}

    if(!bitseqR) return -1;
    if(!bitseqRinv) return -1;

    len = bitseqR->length();

	#ifdef PRINT
	cout<<"\n";
	#endif //PRINT

    fclose(input);
    fclose(hierarchy);
	return ret;
}

int Theorem2::size (){
    return sizeof(Theorem2) + th1->size() + bitseqR->size() + bitseqRinv->size();
}

unsigned int Theorem2::bitsRequired (){
    unsigned int bitsReq = th1->bitsRequired();
    //bits for R and Rinv: (#int) bitmapR + (#int) rank&select overhead (5%) + variable: size
    bitsReq += 2*((bitseqR->length()/W+1 + bitseqR->length()/(W*FACTOR)+1 +1)*W);
    return bitsReq;
}

#endif // THEOREM2_H_INCLUDED
