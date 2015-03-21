/* hutucker.h
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

#ifndef HUTUCKER_H_INCLUDED
#define HUTUCKER_H_INCLUDED

#include <iostream>

#include "permutation.h"
#include "binarytrie.h"

using namespace std;

/** Auxiliar class to build a OABT using a variant of Hu-Tucker algorithm showed in [1].
 *
 *  [1] D. E. Knuth. Art of Computer Programming, Vol. 3 (2nd Edition)
 *
 *  @author Carlos Bedregal
 */

template <class T>
class HuTucker{
    public:
    BNode<T>** seq;
    BNode<T>* root;
    int* levels;
    int start;
    int end;
    int len;
    int weight;

    public:
    HuTucker(Permutation<T>* p);
    HuTucker(T* array, int szArray);
    ~HuTucker();
    BNode<T>* merge(BNode<T>* l, BNode<T>* r);
    int getPosMin();
    void remove(int &pos);
    int findPosCompatible(int& pmin);
    void combination();
    void levelAssignment();
    void recLevelAssign(BNode<T>* curr);
    void recombination();
    void print();
    void printLevels();
};

template <class T>
HuTucker<T>::HuTucker(Permutation<T>* p){
    #ifdef PRINT
        cout<<"- Building HuTucker: "<<p->ro<<" runs\n";
    #endif //PRINT
    start=0;
    root=0;
    weight=0;
    end=p->ro-1;
    len=p->ro;
    seq=new BNode<T>*[p->ro];
    levels=new int[p->ro];

    for(int i=start, point=start; i<=end; i++){
        seq[i]=new BNode<T>(p->Runs[i],i,0,&(p->Runs[i]));
        seq[i]->setEndpoints(point,point+(p->Runs[i])-1);
        point=point+(p->Runs[i]);
    }

    combination();
    levelAssignment();
    recombination();
}

template <class T>
HuTucker<T>::HuTucker(T* array, int szArray){
    #ifdef PRINT
        cout<<"- Building HuTucker: "<<szArray<<" runs ";
    #endif //PRINT
    start=0;
    root=0;
    weight=0;
    end=szArray-1;
    len=szArray;
    seq=new BNode<T>*[szArray];
    levels=new int[szArray];

    for(int i=start, point=start; i<=end; i++){
        seq[i]=new BNode<T>(array[i],i,0,&(array[i]));
        seq[i]->setEndpoints(point,point+(array[i])-1);
        point=point+(array[i]);
    }

    combination();
    levelAssignment();
    recombination();
}

template <class T>
HuTucker<T>::~HuTucker(){
    for(int i=start; i<=end; i++)
        delete seq[i];
    delete[]seq;
    delete[]levels;
}

template <class T>
BNode<T>* HuTucker<T>::merge(BNode<T>* l, BNode<T>* r){
    BNode<T>* n = new BNode<T>(l->w+r->w,l->pos,1,0);
    n->children[0] = l;
    n->children[1] = r;
    n->setEndpoints(l->endpoint[0],r->endpoint[1]);
    return n;
}

template <class T>
int HuTucker<T>::getPosMin(){
    int min=start;
    //sequential search of the minimum element
    for(int i=min+1; i<=end; i++)
        if(seq[min]->w>seq[i]->w)
            min=i ;
    return min;
}

template <class T>
void HuTucker<T>::remove(int &pos){
    seq[pos]=0;
    //sequential update of elements
    for(int i=pos; i<end; i++)
       seq[i]=seq[i+1];
    seq[end]=0;
    end--;
}

template <class T>
int HuTucker<T>::findPosCompatible(int& pmin){
    int pcom, minleft, minright;
    minleft=minright=pmin;
    //find towards left end of pmin's huffman seq
    if(pmin!=start){
        minleft=pmin-1;
        for(int left=pmin-2; left>=start && seq[left+1]->type!=0; left--){
            if(seq[minleft]->w >= seq[left]->w) minleft=left;
        }
    }
    //find towards right end of pmin's huffman seq
    if(pmin!=end){
        minright=pmin+1;
        for(int right=pmin+2; right<=end && seq[right-1]->type!=0; right++){
            if(seq[minright]->w > seq[right]->w) minright=right;
        }
    }
    //min{minleft,minright}
    if(minleft==pmin) pcom=minright;
    else{
        if(minright==pmin) pcom=minleft;
        else pcom=(seq[minleft]->w<=seq[minright]->w) ? minleft : minright;
    }
    return pcom;
}

template <class T>
void HuTucker<T>::combination(){
    int pmin=0, pcom, temp, lastEnd=end;
    while(start!=end){
        pmin = getPosMin();
        pcom = findPosCompatible(pmin);

        if(pmin > pcom){ temp=pmin; pmin=pcom; pcom=temp; }

        seq[pmin]=merge(seq[pmin],seq[pcom]);
        remove(pcom);
    }
    end=lastEnd;
    root=seq[pmin];
}

template <class T>
void HuTucker<T>::levelAssignment(){
    recLevelAssign(root);
}

template <class T>
void HuTucker<T>::recLevelAssign(BNode<T>* curr){
    static int lvl=0;
    if(curr->type==1){//internal node
        weight++; //counter of internal nodes
        lvl++;
        recLevelAssign(curr->children[0]);
        recLevelAssign(curr->children[1]);
        delete curr;
        lvl--;
    }
    else{//leaf node
        levels[curr->pos]=lvl;
        seq[curr->pos]=curr;
    }
}

template <class T>
void HuTucker<T>::recombination(){
    int* stack = new int[len];
    int pos=0;
    int cont=0;
    stack[pos]=cont++;
    while(cont<len){
        pos++;
        stack[pos]=cont++;
        while(pos>0 && levels[stack[pos]]==levels[stack[pos-1]]){//merge
            pos--;
            seq[stack[pos]]=merge(seq[stack[pos]],seq[stack[pos+1]]);
            seq[stack[pos+1]]=0;
            levels[stack[pos]]--;
        }
    }
    root=seq[stack[pos]];
    delete[]stack;
}

template <class T>
void HuTucker<T>::print(){
    int i;
    cout<<"Seq["<<len<<":"<<start<<","<<end<<"]\n";
    for(i=0; i<len; i++)
        if(seq[i]) seq[i]->print();
    cout<<endl;
    if(root!=0){
        root->recPrint();
    }
    cout<<"\n";
}

template <class T>
void HuTucker<T>::printLevels(){
    int i;
    cout<<"Levels: ";
    for(i=0; i<len; i++)
        cout<<levels[i]<<" ";
    cout<<endl;
}


#endif // HUTUCKER_H_INCLUDED
