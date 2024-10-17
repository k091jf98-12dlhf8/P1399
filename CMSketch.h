#ifndef _CMSKETCH_H
#define _CMSKETCH_H

#include <algorithm>
#include <cstring>
#include <string.h>
#include "definition.h"
#include "BOBHash32.h"
#include <iostream>
#include "para.h"

using namespace std;

class CMSketch
{	
private:
	BOBHash32 ** bobhash;
	int* index;
	int**counter;
	int w, d;
	int MAX_CNT;
	int counter_index_size;
	uint64_t hash_value;

public:
	CMSketch(int _w, int _d)
	{
		bobhash = new BOBHash32*[d];
        index = new int[d];
        counter = new int*[d];
		
		counter_index_size = 20;
		w = _w;
		d = _d;
		
		for(int i = 0; i < d; i++)	
		{
			counter[i] = new int[w];
			memset(counter[i], 0, sizeof(int) * w);
		}

		MAX_CNT = (1 << COUNTER_SIZE) - 1;

		for(int i = 0; i < d; i++)
		{
			bobhash[i] = new BOBHash32(i + 1000);
		}
	}
	pair<uint32_t,double> Insert(uint32_t fp,int p)
	{
		int min_cnt = 65535;
		
		for(int i = 0; i < d; i++)
		{
			index[i] = (bobhash[i]->run((const char *)&fp, 4)) % w;
			if(counter[i][index[i]] != MAX_CNT)
			{
				counter[i][index[i]]++;
			}
			if(counter[i][index[i]]<min_cnt)min_cnt=counter[i][index[i]];

		}
		if(p>0){  //a PI flow, check the den
			double den = (double)min_cnt/p;
			if(den<D_thr){return {fp,den};}
		}
		return {fp,-1};
	}
	int Query(uint32_t fp)
	{
		int min_value = MAX_CNT;
		int temp;
		
		for(int i = 0; i < d; i++)
		{
			index[i] = (bobhash[i]->run((const char *)&fp, 4)) % w;
			temp = counter[i][index[i]];
			min_value = temp < min_value ? temp : min_value;
		}
		return min_value;
	
	}
	void Delete(uint32_t fp)
	{
		for(int i = 0; i < d; i++)
		{
			index[i] = (bobhash[i]->run((const char *)&fp, 4)) % w;
			counter[i][index[i]] --;
		}
	}
	~CMSketch()
	{
		for(int i = 0; i < d; i++)	
		{
			delete []counter[i];
		}


		for(int i = 0; i < d; i++)
		{
			delete bobhash[i];
		}
	}
};
#endif//_CMSKETCH_H