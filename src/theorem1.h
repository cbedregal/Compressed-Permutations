/* theorem1.h
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

#ifndef THEOREM1_H_INCLUDED
#define THEOREM1_H_INCLUDED

#include"theorem.h"

/** Implementation of Compressed Data Structure for Permutations based on 
 *	Wavelet Tree (Hu-Tucker) using ascending sub-sequences (Runs).
 *  (Theorem2 or TH2 for practical use)
 *
 *  [1] J. Barbay and G. Navarro, Compressed Representation of Permutations,
 *  and Applications.
 *
 *  @author Carlos Bedregal
 */

class Theorem1:public Theorem{
    public:
    WaveletTree<int> *wt;
    //int waste;

    public:
    Theorem1();
    Theorem1(Permutation<int> *p);
    virtual ~Theorem1();

    WaveletTree<int> * tree();

    uint pi(int i);
    uint recPi(WTNode* node, WTNode* parent, int j);
    uint piInv(int i);
    uint recPiInv(WTNode* node, int i);

    int save (char* fname);
    int recSave(WTNode* node, FILE* fp, uint* shape, uint& curr);

    int load (char* fname);
    int loadWT(FILE* fp,FILE* fh);
    int recLoad(FILE* fp, WTNode* node, uint* shape, uint& curr);

    int size();
    void recSize(WTNode* node, int& size);

    unsigned int bitsRequired();
};

Theorem1::Theorem1(){
    wt=0;
}

Theorem1::Theorem1(Permutation<int> *p){
    assert(p!=0);
    assert(p->len>0);
    #ifdef PRINT
        cout<<"+ Construyendo Th1 en "<<p->len<<" elementos\n";
    #endif //PRINT
    //Permutation<int>* p = new Permutation<int>(array,n);

    len=p->len;
    wt=new WaveletTree<int>(p->array,p->Runs,p->ro);
    cout<<"nodes: "<<wt->weight<<endl;
}

Theorem1::~Theorem1(){
    delete wt;
}

WaveletTree<int>* Theorem1::tree(){
    return wt;
}

uint Theorem1::pi(int i){
    return recPi(wt->root,0,i+1);
}

uint Theorem1::recPi(WTNode* node, WTNode* parent, int j){
    static int s;

    //s=node->size;
    s=node->bitseq->length();

    #ifdef DEBUG
        cout<<"\tDOWN: nodo: "<<node->size<<", s: "<<s<<", j: "<<j<<", rank0(B,s-1): "<<node->bitseq->rank0(s-1)<<endl;
    #endif //DEBUG

    //downward traversal to determine leaf v and offset j
    //a) go down to the left
    if(node->bitseq->rank0(s-1) >= (unsigned int)j){
        if(!node->children[0])
            j=node->bitseq->select0(j)+1;
        else
            j=recPi(node->children[0],node,j);
    }
    //b) go down to the right
    else{
        j=j-node->bitseq->rank0(s-1);
        if(!node->children[1])
            j=node->bitseq->select1(j)+1;
        else
            j=recPi(node->children[1],node,j);
    }

    #ifdef DEBUG
        cout<<"\tUP: nodo: "<<node->size<<", s: "<<s<<", j: "<<j;
    #endif //DEBUG

    //we've reach the root
    if(!parent)
        return j-1;

    #ifdef DEBUG
        cout<<", left_child?: "<<(node==parent->children[0])<<endl;
    #endif //DEBUG

    //upward traversal of nodes in the recursion stack
    //a) left child of parent
    if(node==parent->children[0])
        j=parent->bitseq->select0(j);
    //b) right child of parent
    else
        j=parent->bitseq->select1(j);

    return ++j;
}

uint Theorem1::piInv(int i){
    return recPiInv(wt->root,i);
}

uint Theorem1::recPiInv(WTNode* node, int i){
    static int p=0;
    static int s;

    //is leaf?
    if(!node){
        int result=p+i;
        p=0;
        return result;
    }

    //s=node->size;
    s=node->bitseq->length();

    #ifdef DEBUG
        cout<<"\tnodo: "<<node->size<<", s: "<<s<<", i: "<<i<<", p: "<<p<<", B[i]: "<<node->bitseq->access(i)<<endl;
    #endif //DEBUG

    //B[i]=1, go down to the right
    if(node->bitseq->access(i)){
        p=p+node->bitseq->rank0(s-1);
        return recPiInv(node->children[1],node->bitseq->rank1(i)-1);
    }
    //B[i]=0, go down to the left
    else{
        return recPiInv(node->children[0],node->bitseq->rank0(i)-1);
    }
}

/* saves structure TH1's bitmaps into files with prefix "fname" through recursive
 * method recSave
 * - fname: stores th1's bitsequences
 * - fname.idx: stores the tree shape
 */
int Theorem1::save (char* fname){
	char fname2[128];
	strcpy(fname2,fname);
	strcat(fname2,".idx");
	uint* shape = new uint[uint_len(wt->weight*2,1)+1];
	uint curr=0;

    FILE * output;
    output = fopen(fname,"wb");

	//save each node recursively
    int ret = recSave(wt->root,output,shape,curr);
    fclose(output);

	FILE * hierarchy;
	hierarchy = fopen(fname2,"wb");

	//save tree structure
	fwrite(&curr,sizeof(uint),1,hierarchy);
	fwrite(shape,sizeof(uint),uint_len(wt->weight*2,1)+1,hierarchy);
	fclose(hierarchy);

	delete[]shape;
    return ret;
}

int Theorem1::recSave(WTNode* node, FILE* fp, uint* shape, uint& curr){
    //static char child;
	static bool exit;
    if(!node){
        //child='0';
        //fwrite(&child,sizeof(char),1,fp);
		bitclean(shape,curr++);
        return 0;
    }
    else{
        //child='1';
        //fwrite(&child,sizeof(char),1,fp);
		bitset(shape,curr++);
		//save node's bitsequence+headers into fp
		exit = node->save(fp);
        assert(exit==0);
        return recSave(node->children[0],fp,shape,curr) + recSave(node->children[1],fp,shape,curr);
    }
}

/* loads in memory a previously saved TH1 structure from file "fname" through
 * method loadWT.
 * - fname contains the bitsequence of each node
 * - fname.idx: contains the three shape
 */
int Theorem1::load(char* fname){
	char fname2[128];
	strcpy(fname2,fname);
	strcat(fname2,".idx");

    FILE * input;
    input = fopen(fname,"rb");
	FILE * hierarchy;
	hierarchy = fopen(fname2,"rb");

    if(loadWT(input,hierarchy)==-1){
		cout<<"@Theorem1::load()\n";
		return -1;
	}

	#ifdef PRINT
	cout<<"\n";
	#endif //PRINT

    fclose(input);
	fclose(hierarchy);
    return 0;
}

/* loads in memory a previously saved TH1 through recursive method recLoad.
 * - fp contains the bitsequence of each node
 * - fh contains the three shape
 */
int Theorem1::loadWT(FILE* fp, FILE* fh){
	//load tree structure fro fh
	uint curr=0;
	if(fread(&curr,sizeof(uint),1,fh)!=1 || curr==0){
		cout<<"@Theorem1::loadWT(): fread(tree.size)\n";
		return -1;
	}

	uint* shape = new uint[uint_len(curr,1)];
	if (fread (shape,sizeof(uint),uint_len(curr,1),fh) != uint_len(curr,1)){
		cout<<"@Theorem1::loadWT(): fread(tree.shape)\n";
		return -1;
	}

	//load root of wavelet tree
	curr=0;
	wt = new WaveletTree<int>();
    wt->root = new WTNode(); wt->weight++;
    if(bitget(shape,curr++)!=1){
		cout<<"@Theorem1::loadWT(): root flag\n";
		return -1;
	}
    if(wt->root->load(fp)!=0){
		cout<<"@Theorem1::loadWT() root->load\n";
		return -1;
	}
	len = wt->root->bitseq->length();

	//load wavelet tree recursively
    return recLoad(fp,wt->root,shape,curr);
}

int Theorem1::recLoad(FILE* fp, WTNode* node, uint* shape, uint& curr){
    static bool child;
    //left child
	child = bitget(shape,curr++);
    if(child){ //next to read is child of node
        node->children[0] = new WTNode(); wt->weight++;
        if(node->children[0]->load(fp)!=0){
			cout<<"@Theorem1::recLoad() left_child->load\n";
			return -1;
		}
        recLoad(fp,node->children[0],shape,curr);
    }
    else{ //no child
		node->children[0] = 0;
    }
    //right child
    child = bitget(shape,curr++);
    if(child){ //next to read is child of node
        node->children[1] = new WTNode(); wt->weight++;
        if(node->children[1]->load(fp)!=0){
			cout<<"@Theorem1::recLoad() right_child->load\n";
			return -1;
		}
		recLoad(fp,node->children[1],shape,curr);
    }
    else{ //no child
		node->children[1] = 0;
    }

    return 0;
}

int Theorem1::size (){
    int size = 0;
    //waste = 0;
    recSize(wt->root,size);
    return sizeof(Theorem1) + sizeof(WaveletTree<int>) + size;
}

void Theorem1::recSize(WTNode* node, int& size){
    size += node->size();
    //waste += uint_len(node->bitseq->length(),1)*32-node->bitseq->length();
    if(node->children[0]) recSize(node->children[0],size);
    if(node->children[1]) recSize(node->children[1],size);
}

unsigned int Theorem1::bitsRequired(){
    unsigned int b=0;
    wt->recBitsRequired(wt->root,b);
	//cout<<"#B: "<<b<<endl;
	assert(b>0);
    return b;
}

#endif // THEOREM1_H_INCLUDED
