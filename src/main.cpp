#include<iostream>
#include<fstream>
#include<vector>

#define PRINT

#define BRW 0
#define RRRL 1
#define RRR 2

int bitseqFlag=BRW;

#include"theorem1.h"
#include"theorem2.h"

using namespace std;

void printArray(int* array, int n){
    cout<<"["<<n<<"]: ";
    for(int i=0; i<n;  i++){
        cout<<array[i]<<" ";
    }
    cout<<endl;
}

int* createArray(int ro, int* runs,int &size){
    int *array=0;
    int i,j,k,last;

    if(runs==0){//random
        ro = rand()%11+10;
        runs=new int[ro];
        for(i=0; i<ro; i++)
            runs[i]=rand()%10+1;
    }

    for(i=0,size=0; i<ro; size+=runs[i++]);
    array=new int[size];

    for(i=0,k=0,last=size; i<ro;  i++){
        last=last-runs[i];
        for(j=0; j<runs[i]; j++,k++)
            array[k]=last+j;
    }

    return array;
}

int main(int argc, char* argv[]){

	int *array, size;
    int runs[]={5,2,7,2,1,1,1,2,4,5};

    int ro=sizeof(runs)/sizeof(int);
    array=createArray(ro,runs,size);
    printArray(array,size);

    Permutation <int> p (array,size);
    p.findRuns();

    Theorem *th = new Theorem1(&p);
    //Theorem *th = new Theorem2(&p);

	cout<<th->pi(0)<<endl;

	//th->save("test.x");
	//th->load("test.x");

	cout<<th->piInv(0)<<endl;

    return 0;
}
