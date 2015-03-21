/* permutation.h
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

#ifndef PERMUTATION_H_INCLUDED
#define PERMUTATION_H_INCLUDED

#include<iostream>

using namespace std;

/** Auxiliar class to handle permutations and the identification of runs within it.
 *
 *  @author Carlos Bedregal
 */

template <class T>
class Permutation{
    public:
    T* array;
    unsigned int len;
    int ro;
    int* Runs;
    int tau;
    int* SRuns;
    int Hro;
    int* HRuns;

    public:
    Permutation(unsigned int n);
    Permutation(T* A, unsigned int n);
    ~Permutation();
    T* operator[] (int pos);
    void print();
    void findRuns();
    void findSRuns();
    void findHRuns();
	void copyArray(T* perm);
};

template <class T>
Permutation<T>::Permutation(unsigned int n){
    len=n;
    ro=0;
    tau=0;
    Hro=0;
    //Runs=new int[n];
    //for(int i=0;i<n;Runs[i]=0,i++);
}

template <class T>
Permutation<T>::Permutation(T* A, unsigned int n){
    #ifdef PRINT
        cout<<"- Building Permutation: "<<n<<" elementos\n";
    #endif //PRINT

    len=n;
    copyArray(A);
    ro=0;
    tau=0;
    Hro=0;
}

template <class T>
Permutation<T>::~Permutation(){
    if(ro!=0) delete[]Runs;
    if(tau!=0) delete[]SRuns;
    if(Hro!=0) delete[]HRuns;
}

template <class T>
T* Permutation<T>::operator[] (int pos){
    return &(array[pos]);
}

template <class T>
void Permutation<T>::copyArray(T* perm){
    array=perm;
}

template <class T>
void Permutation<T>::print(){
    int i;
    cout<<"Array["<<len<<"]\n";
    for(i=0; i<len; i++)
        cout<<array[i]<<" ";
    if(ro!=0){
        cout<<"\nRuns["<<ro<<"]\n";
        for(i=0; i<ro; i++)
            cout<<Runs[i]<<" ";
    }
    if(tau!=0){
        cout<<"\nSRuns["<<tau<<"]\n";
        for(i=0; i<tau; i++)
            cout<<SRuns[i]<<" ";
    }
    if(Hro!=0){
        cout<<"\nHRuns["<<Hro<<"]\n";
        for(i=0; i<Hro; i++)
            cout<<HRuns[i]<<" ";
    }
    cout<<"\n\n";
}

template <class T>
void Permutation<T>::findRuns(){
	if(ro!=0) return;
	#ifdef PRINT
		cout<<"\t- Identifying Runs\n";
	#endif //PRINT
    unsigned int i;
    for(i=1; i<len; i++){
        if(array[i]<array[i-1]) //run down_step
            ro++;
    }
    ro++;

    Runs=new int[ro];
    for(int j=0;j<ro;Runs[j]=0,j++);

    ro=0;
    for(i=1; i<len; i++){//calculate the size of each run
        Runs[ro]++;
        if(array[i]<array[i-1]) //run down_step
            ro++;
    }
    Runs[ro++]++; //last position
}

template <class T>
void Permutation<T>::findSRuns(){
	if(tau!=0) return;
	#ifdef PRINT
		cout<<"\t- Identifying SRuns\n";
	#endif //PRINT
    unsigned int i;
    for(i=1; i<len; i++){
        if(array[i]!=array[i-1]+1)//strict run down_step
            tau++;
    }
	tau++;

    SRuns=new int[tau];
    for(int j=0;j<tau;SRuns[j]=0,j++);

    tau=0;
    for(i=1; i<len; i++){//calculate the size of each run
        SRuns[tau]++;
        if(array[i]!=array[i-1]+1)//strict run down_step
            tau++;
    }
    SRuns[tau++]++; //last position
}

template <class T>
void Permutation<T>::findHRuns(){
	if(Hro!=0) return;
	#ifdef PRINT
		cout<<"\t- Identifying HRuns\n";
	#endif //PRINT
    int i;
    for(i=1, Hro=0; i<tau; i++){
        if(SRuns[i]<SRuns[i-1]) //SRun down_step
            Hro++;
    }
    Hro++;

    HRuns=new int[Hro];
    for(i=0;i<Hro;HRuns[i]=0,i++);

    Hro=0;
    for(i=1; i<tau; i++){//calculate the size of each SRun
        HRuns[Hro]++;
        if(SRuns[i]<SRuns[i-1]) //SRun down_step
            Hro++;
    }
    HRuns[Hro++]++; //last position
}
#endif // PERMUTATION_H_INCLUDED
