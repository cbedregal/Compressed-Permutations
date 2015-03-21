/* static_bitsequence_brw32.cpp
   Copyright (C) 2005, Rodrigo Gonzalez, all rights reserved.

   New RANK, SELECT, SELECT-NEXT and SPARSE RANK implementations.

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

#include "static_bitsequence_brw32.h"
#include <cassert>
#include <cmath>
// #include <sys/types.h>


/////////////
//Rank(B,i)//
/////////////
//_factor = 0  => s=W*lgn
//_factor = P  => s=W*P
//Is interesting to notice
//factor=2 => overhead 50%
//factor=3 => overhead 33%
//factor=4 => overhead 25%
//factor=20=> overhead 5%

static_bitsequence_brw32::static_bitsequence_brw32(){
  data=NULL;
//  this->owner = true;
  this->len=0;
  //this->factor=0;
}

static_bitsequence_brw32::static_bitsequence_brw32( uint *bitarray, uint _n, uint _factor){
  /*cout << "*****" << endl;
  cout << bitarray << endl;
  cout << _n << endl;
  cout << _factor << endl; */
	this->len = _n;
	this->ones = len/W+1;
  if(_factor==0) exit(-1);
  data=new uint[_n/W+1];
  for(uint i=0;i<uint_len(_n,1);i++)
    data[i] = bitarray[i];
  for(uint i=uint_len(_n,1);i<_n/W+1;i++)
    data[i] = 0;
  //this->owner = true;
  //this->n=_n;
  //uint lgn=bits(len-1);
  //this->factor=_factor;
  //if (_factor==0) this->factor=lgn;
  //else this->factor=_factor;
  //b=32;
  //s=b*this->factor;
  //integers = n/W+1;
  BuildRank();
  
}

static_bitsequence_brw32::~static_bitsequence_brw32() {
  delete [] Rs;
  delete [] data;
}

//Metodo que realiza la busqueda d
void static_bitsequence_brw32::BuildRank(){
  uint num_sblock = len/S;
  Rs = new uint[num_sblock+5];// +1 pues sumo la pos cero
  for(uint i=0;i<num_sblock+5;i++)
    Rs[i]=0;
  uint j;
  Rs[0]=0;
  for (j=1;j<=num_sblock;j++) {
    Rs[j]=Rs[j-1];
    Rs[j]+=BuildRankSub((j-1)*FACTOR,FACTOR);
  }
}

uint static_bitsequence_brw32::BuildRankSub(uint ini,uint bloques){
  uint rank=0,aux;
  for(uint i=ini;i<ini+bloques;i++) {
    if (i < this->ones) {
      aux=data[i];
      rank+=popcount(aux);
    }
  }
  return rank; //retorna el numero de 1's del intervalo

}

uint static_bitsequence_brw32::rank1(uint i) {
  ++i;
  uint resp=Rs[i/S];
  uint aux=(i/S)*FACTOR;
  for (uint a=aux;a<i/W;a++)
    resp+=popcount(data[a]);
  resp+=popcount(data[i/W]  & ((1<<(i & mask31))-1));
  return resp;
}

bool static_bitsequence_brw32::access(uint i) {
  return (1u << (i % W)) & data[i/W];
}

int static_bitsequence_brw32::save(FILE *f) {
  uint wr = BRW32_HDR;
  if (f == NULL) return 20;
  if( fwrite(&wr,sizeof(uint),1,f) != 1 ) return -1;
  if (fwrite (&len,sizeof(uint),1,f) != 1) return 21;
  if (fwrite (&wr,sizeof(uint),1,f) != 1) return 21;
  if (fwrite (data,sizeof(uint),len/W+1,f) != len/W+1) return 21;
  if (fwrite (Rs,sizeof(uint),len/S+1,f) != len/S+1) return 21;
  return 0;
}

static_bitsequence_brw32 * static_bitsequence_brw32::load(FILE *f) {
  if (f == NULL) return NULL;
  uint type;
  if(fread(&type,sizeof(uint),1,f)!=1) return NULL;
  if(type!=BRW32_HDR) { cout << "type:"<<type<<endl; return NULL; }
  static_bitsequence_brw32 * ret = new static_bitsequence_brw32();
  if (fread (&ret->len,sizeof(uint),1,f) != 1) return NULL;
  //ret->b=32; // b is a word
  if (fread (&type,sizeof(uint),1,f) != 1) return NULL;
  //ret->s=ret->b*ret->factor;
  uint aux=(ret->len+1)%W;
  if (aux != 0)
    ret->ones = (ret->len+1)/W+1;
  else
    ret->ones = (ret->len+1)/W;
  ret->data = new uint[ret->len/W+1];
  if (!ret->data) return NULL;
  if (fread (ret->data,sizeof(uint),ret->len/W+1,f) != ret->len/W+1) return NULL;
  ret->Rs= new uint[ret->len/S+1];
  if (!ret->Rs) return NULL;
  if (fread (ret->Rs,sizeof(uint),ret->len/S+1,f) != ret->len/S+1) return NULL;
	ret->len = ret->len;
  return ret;
}

uint static_bitsequence_brw32::SpaceRequirementInBits() {
  return uint_len(len,1)*sizeof(uint)*8+(len/S)*sizeof(uint)*8;
}

uint static_bitsequence_brw32::size() {
  return sizeof(static_bitsequence_brw32)+SpaceRequirementInBits()/8;
}

uint static_bitsequence_brw32::SpaceRequirement() {
  return len/8+(len/S)*sizeof(uint)+sizeof(static_bitsequence_brw32);
}

uint static_bitsequence_brw32::prev2(uint start) {
      // returns the position of the previous 1 bit before and including start.
      // tuned to 32 bit machine

      uint i = start >> 5;
      int offset = (start % W);
      uint answer = start;
      uint val = data[i] << (Wminusone-offset);

      if (!val) { val = data[--i]; answer -= 1+offset; }

      while (!val) { val = data[--i]; answer -= W; }

      if (!(val & 0xFFFF0000)) { val <<= 16; answer -= 16; }
      if (!(val & 0xFF000000)) { val <<= 8; answer -= 8; }

      while (!(val & 0x80000000)) { val <<= 1; answer--; }
      return answer;
}

uint static_bitsequence_brw32::prev(uint start) {
      // returns the position of the previous 1 bit before and including start.
      // tuned to 32 bit machine

      uint i = start >> 5;
      int offset = (start % W);
      uint aux2 = data[i] & (-1u >> (31-offset));

      if (aux2 > 0) {
                if ((aux2&0xFF000000) > 0) return i*W+23+prev_tab[(aux2>>24)&0xFF];
                else if ((aux2&0xFF0000) > 0) return i*W+15+prev_tab[(aux2>>16)&0xFF];
                else if ((aux2&0xFF00) > 0) return i*W+7+prev_tab[(aux2>>8)&0xFF];
                else  return i*W+prev_tab[aux2&0xFF]-1;
      }
      for (uint k=i-1;;k--) {
         aux2=data[k];
         if (aux2 > 0) {
                if ((aux2&0xFF000000) > 0) return k*W+23+prev_tab[(aux2>>24)&0xFF];
                else if ((aux2&0xFF0000) > 0) return k*W+15+prev_tab[(aux2>>16)&0xFF];
                else if ((aux2&0xFF00) > 0) return k*W+7+prev_tab[(aux2>>8)&0xFF];
                else  return k*W+prev_tab[aux2&0xFF]-1;
         }
      }
      return 0;
}

uint static_bitsequence_brw32::next(uint k) {
        uint count = k;
        uint des,aux2;
        des=count%W;
        aux2= data[count/W] >> des;
        if (aux2 > 0) {
                if ((aux2&0xff) > 0) return count+select_tab[aux2&0xff]-1;
                else if ((aux2&0xff00) > 0) return count+8+select_tab[(aux2>>8)&0xff]-1;
                else if ((aux2&0xff0000) > 0) return count+16+select_tab[(aux2>>16)&0xff]-1;
                else {return count+24+select_tab[(aux2>>24)&0xff]-1;}
        }

        for (uint i=count/W+1;i<this->ones;i++) {
                aux2=data[i];
                if (aux2 > 0) {
                        if ((aux2&0xff) > 0) return i*W+select_tab[aux2&0xff]-1;
                        else if ((aux2&0xff00) > 0) return i*W+8+select_tab[(aux2>>8)&0xff]-1;
                        else if ((aux2&0xff0000) > 0) return i*W+16+select_tab[(aux2>>16)&0xff]-1;
                        else {return i*W+24+select_tab[(aux2>>24)&0xff]-1;}
                }
        }
        return len;
}

uint static_bitsequence_brw32::select1(uint x) {
  // returns i such that x=rank(i) && rank(i-1)<x or n if that i not exist
  // first binary search over first level rank structure
  // then sequential search using popcount over a int
  // then sequential search using popcount over a char
  // then sequential search bit a bit

  //binary search over first level rank structure
  uint l=0, r=len/S;
  uint mid=(l+r)/2;
  uint rankmid = Rs[mid];
  while (l<=r) {
    if (rankmid<x)
      l = mid+1;
    else
      r = mid-1;
    mid = (l+r)/2;
    rankmid = Rs[mid];
  }
  //sequential search using popcount over a int
  uint left;
  left=mid*FACTOR;
  x-=rankmid;
        uint j=data[left];
        uint ones = popcount(j);
        while (ones < x) {
    x-=ones;left++;
    if (left > this->ones) return len;
          j = data[left];
      ones = popcount(j);
        }
  //sequential search using popcount over a char
  left=left*B;
  rankmid = popcount8(j);
  if (rankmid < x) {
    j=j>>8;
    x-=rankmid;
    left+=8;
    rankmid = popcount8(j);
    if (rankmid < x) {
      j=j>>8;
      x-=rankmid;
      left+=8;
      rankmid = popcount8(j);
      if (rankmid < x) {
        j=j>>8;
        x-=rankmid;
        left+=8;
      }
    }
  }

  // then sequential search bit a bit
        while (x>0) {
    if  (j&1) x--;
    j=j>>1;
    left++;
  }
  return left-1;
}

uint static_bitsequence_brw32::select0(uint x) {
  // returns i such that x=rank_0(i) && rank_0(i-1)<x or n if that i not exist
  // first binary search over first level rank structure
  // then sequential search using popcount over a int
  // then sequential search using popcount over a char
  // then sequential search bit a bit

  //binary search over first level rank structure
  if(x==0) return 0;
  uint l=0, r=len/S;
  uint mid=(l+r)/2;
  uint rankmid = mid*FACTOR*W-Rs[mid];
  while (l<=r) {
    if (rankmid<x)
      l = mid+1;
    else
      r = mid-1;
    mid = (l+r)/2;
    rankmid = mid*FACTOR*W-Rs[mid];
  }
  //sequential search using popcount over a int
  uint left;
  left=mid*FACTOR;
  x-=rankmid;
  uint j=data[left];
  uint zeros = W-popcount(j);
  while (zeros < x) {
    x-=zeros;left++;
    if (left > this->ones) return len;
    j = data[left];
    zeros = W-popcount(j);
  }
  //sequential search using popcount over a char
  left=left*B;
  rankmid = 8-popcount8(j);
  if (rankmid < x) {
    j=j>>8;
    x-=rankmid;
    left+=8;
    rankmid = 8-popcount8(j);
    if (rankmid < x) {
      j=j>>8;
      x-=rankmid;
      left+=8;
      rankmid = 8-popcount8(j);
      if (rankmid < x) {
        j=j>>8;
        x-=rankmid;
        left+=8;
      }
    }
  }

  // then sequential search bit a bit
  while (x>0) {
    if  (j%2 == 0 ) x--;
    j=j>>1;
    left++;
  }
  left--;
  if (left > len)  return len;
  else return left;
}
