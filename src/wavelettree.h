/* wavelettree.h
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

#ifndef WAVELETTREE_H_INCLUDED
#define WAVELETTREE_H_INCLUDED

#include "hutucker.h"
#include "waveletnode.h"

using namespace std;

/** Class for wavelet tree data structure. Builds a wavelet tree form a Hu-Tucker shaped binarytrie,
 *  it also sorts (merging the nodes) the original permutation.
 *
 *  @author Carlos Bedregal
 */

template <class T>
class WaveletTree{
    public:
    T* array; //array of values
    //HuTucker<T>* ht; //pointer to HuTucker tree
    WTNode* root;
    int weight;

    public:
    WaveletTree();
    WaveletTree(HuTucker<T>* ht, T* array);
    WaveletTree(T* array, int* runs, int ro);
    ~WaveletTree();
    void recBuild(BNode<T>* bNode, WTNode* wNode);
    void recPrint(WTNode* node);
    void recDestruct(WTNode* node);
    void recBitsRequired(WTNode* node, unsigned int& bitsReq);
};

template <class T>
WaveletTree<T>::WaveletTree(){
    array=0;
    root=0;
    weight=0;
}

template <class T>
WaveletTree<T>::WaveletTree(HuTucker<T>* ht, T* array){
    assert(ht!=0);
    assert(array!=0);
    this->array=array;
    root=new WTNode(ht->root->w);
    weight=1;
    recBuild(ht->root, root);
    assert(weight==ht->weight);
}

template <class T>
WaveletTree<T>::WaveletTree(T* array, int* runs, int ro){
    #ifdef PRINT
        cout<<"- Building WaveletTree: "<<ro<<" run\n";
    #endif //PRINT
    assert(array!=0);
    assert(runs!=0);
    assert(ro>0);
    this->array=array;

    HuTucker<T>* ht = new HuTucker<T>(runs,ro);
    //HuTucker<T>* ht = new HuTucker<T>(p);
    assert(ht->len=ro);
    assert(ht->root!=0);

    root=new WTNode(ht->root->w);
    weight=1;
    recBuild(ht->root, root);

    assert(weight==ht->weight);
    #ifdef PRINT
        cout<<"- "<<weight<<" nodes\n";
    #endif //PRINT

    delete ht;
}

template <class T>
WaveletTree<T>::~WaveletTree(){
    recDestruct(root);
}

template <class T>
void WaveletTree<T>::recBuild(BNode<T>* bNode, WTNode* wNode){
    //build nodes on the wavelet-tree only for internal nodes of hu-tucker
    if(bNode->children[0]->type!=0){
        weight++;
        wNode->children[0]=new WTNode(bNode->children[0]->w);
        recBuild(bNode->children[0],wNode->children[0]);
    }
    if(bNode->children[1]->type!=0){
        weight++;
        wNode->children[1]=new WTNode(bNode->children[1]->w);
        recBuild(bNode->children[1],wNode->children[1]);
    }

    //merge of nodes (merge both bitmaps and sort the area covered by the nodes)
    int i,j,k;
    uint* bitmap = new uint[uint_len(bNode->w,1)];
    //int* mergeArea=new int[wNode->size];
    int* mergeArea=new int[bNode->w];
    for(i=0,j=bNode->children[1]->w-1,k=0; i<bNode->children[0]->w && j>=0; k++){
        if(array[bNode->endpoint[0]+i]<array[bNode->endpoint[1]-j]){
            //bitclean(wNode->bitmap,k);
            bitclean(bitmap,k);
            mergeArea[k]=array[bNode->endpoint[0]+i];
            i++;
        }
        else{
            //bitset(wNode->bitmap,k);
            bitset(bitmap,k);
            mergeArea[k]=array[bNode->endpoint[1]-j];
            j--;
        }
    }
    for(; i<bNode->children[0]->w; i++,k++){
        //bitclean(wNode->bitmap,k);
        bitclean(bitmap,k);
        mergeArea[k]=array[bNode->endpoint[0]+i];
    }
    for(; j>=0; j--,k++){
        //bitset(wNode->bitmap,k);
        bitset(bitmap,k);
        mergeArea[k]=array[bNode->endpoint[1]-j];
    }
    //for(k=0; k<wNode->size; k++)
    for(k=0; k<bNode->w; k++)
        array[k+bNode->endpoint[0]]=mergeArea[k];
    //end of merge

    //wNode->createBitseq();
    wNode->createBitseq(bitmap,bNode->w);

    #ifdef DEBUG
    //for(k=0; k<root->size; k++)
    //    cout<<array[k]<<" ";
    //cout<<endl;
    #endif //DEBUG

    delete[]bitmap;
    delete[]mergeArea;
}

template <class T>
void WaveletTree<T>::recPrint(WTNode* node){
    if(!node){
        cout<<"null\n"; return;
    }
    node->print();
    cout<<"down\n";
    cout<<"left: "; recPrint(node->children[0]);
    cout<<"right: "; recPrint(node->children[1]);
    cout<<"up\n";
}

template <class T>
void WaveletTree<T>::recDestruct(WTNode* node){
    if(!node) return;
    recDestruct(node->children[0]);
    recDestruct(node->children[1]);
    delete node;
}

template <class T>
void WaveletTree<T>::recBitsRequired(WTNode* node, unsigned int& bitsReq){
	int nodeBits = node->bitseq->length();
	//bits requiered for: [  (#ints) node's bitmap   +   (#ints) rank&select (5%)  ] * W=32
    bitsReq += (nodeBits/W+1 + nodeBits/(W*FACTOR)+1)*W;  //((static_bitsequence_brw32*)node->bitseq)->SpaceRequirementInBits();
    //bits requiered for rank&select constants: size of bitmap         //WARNING: the real implementation considers 3 constants: size+factor+flag
    //bitsReq += 1*W;
    //bits for the tree pointers: 2 bits (one for children: 1 or 0)
    //bitsReq += 2;
    if(node->children[0]) recBitsRequired(node->children[0],bitsReq);
    if(node->children[1]) recBitsRequired(node->children[1],bitsReq);
}

#endif // WAVELETTREE_H_INCLUDED
