/* binarytrie.h
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

#ifndef BINARYTRIE_H_INCLUDED
#define BINARYTRIE_H_INCLUDED

#include<iostream>

using namespace std;

/** Auxiliar class to handle binary nodes.
 *
 *  @author Carlos Bedregal
 */

template <class T>
class BNode{
    public:
    int w; //weight
    int pos; //position
    bool type; //0=external, 1=internal
    int endpoint[2];
    BNode* children[2];
    T* obj;

    BNode();
    BNode(int v, int p, int t, T* o);
    ~BNode(){};
    void print();
    void setEndpoints(int start, int end);
    void recPrint(int level=0);
};

template <class T>
BNode<T>::BNode(){
    w=-1;
    obj=0;
    pos=-1;
    type=0;
    endpoint[0]=endpoint[1]=-1;
    children[0]=children[1]=0;
}

template <class T>
BNode<T>::BNode(int v, int p, int t, T* o){
    w=v;
    obj=o;
    pos=p;
    type=t;
    endpoint[0]=endpoint[1]=-1;
    children[0]=children[1]=0;
}

template <class T>
void BNode<T>::setEndpoints(int start, int end){
    endpoint[0]=start;
    endpoint[1]=end;
}

template <class T>
void BNode<T>::print(){
    cout<<"@:"<<this<<": w="<<w<<" pos="<<pos<<" type="<<type;
    if(obj) cout<<" obj: "<<*obj;
    cout<<" ["<<endpoint[0]<<","<<endpoint[1]<<"]";
    cout<<endl;
}

template <class T>
void BNode<T>::recPrint(int level){
    cout<<"@L:"<<level<<": w="<<w<<" pos="<<pos<<" type="<<type;
    if(obj) cout<<" obj: "<<*obj;
    cout<<" ["<<endpoint[0]<<","<<endpoint[1]<<"]";
    cout<<endl;
    if(children[0]){
        children[0]->recPrint(level+1);
        children[1]->recPrint(level+1);
        cout<<"--------\n";
    }
}


#endif // BINARYTRIE_H_INCLUDED
