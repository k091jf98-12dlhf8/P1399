
#ifndef CLASS2_H
#define CLASS2_H

#include "hash.h"
#include "BOBHash32.h"
#include <cstring>
#include <algorithm>
#include <vector>
#include <stdint.h>
#include "para.h"
#include <time.h>
#include <stdlib.h>
#include <immintrin.h>

using namespace std;

class Bucket_S
{
public:
    S_ID_type *FP;  
    S_F_type *f;  
	S_P_type *p;

    Bucket_S() {
        FP = new S_ID_type[S_entry_num];
        f = new S_F_type[S_entry_num];
		p = new S_P_type[S_entry_num];
        memset(FP, 0, sizeof(S_ID_type) * S_entry_num);
		memset(f, 0, sizeof(S_F_type) * S_entry_num);
		memset(p, 0, sizeof(S_P_type) * S_entry_num);
    }
};



class Bucket_L
{
public:
    L_ID_type *ID;
    L_of_type *of;  
    Bucket_L() {
        ID = new L_ID_type[L_record_num];
        of = new L_of_type[L_record_num];
        memset(ID, 0, sizeof(L_ID_type) * L_record_num);
		memset(of, 0, sizeof(L_of_type) * L_record_num);
    }
};




class L_virtual
{
public:
	virtual int report(uint32_t ID, uint8_t f,uint8_t p,char mode) = 0; 
};




class S_Sketch{
public:
    Bucket_S *Filter;  
    BOBHash32 *P_Hash;
	int w_cnt;

	L_virtual* L_vir;

    S_Sketch() {
        Filter = new Bucket_S[S_bucket_num];
        for (int i = 0; i < S_bucket_num; i++){
            Filter[i] = Bucket_S();
        }
        
        P_Hash = new BOBHash32(97);
		w_cnt = 0;
    }


bool insert(uint32_t key, int count = 1){
    	w_cnt++;
		if (w_cnt == WindowSize) {
			w_cnt = 0;
			__m256i mask = _mm256_set1_epi8(0xFD);	
			for (int i = 0; i < S_bucket_num; i++) {
				__m256i data = _mm256_loadu_si256((__m256i*)Filter[i].p);  
				data = _mm256_and_si256(data, mask);  
				_mm256_storeu_si256((__m256i*)Filter[i].p, data); 
			}
		}
	
		
		S_ID_type keyfp = key%pow(2,S_ID_len);
		int hash_index = P_Hash->run((const char *)&key, 4) % S_bucket_num;
		Bucket_S Bid = Filter[hash_index];

		int match_index = -1; 
		int first_empty=-1;  
		int p_min_pos=-1;  
		int p_min_value=99;   
		
	 
		__m256i target_FP = _mm256_set1_epi32(key); 
		__m256i target_empty = _mm256_set1_epi32(0);
   		for (int i = 0; i < 32; i += 8) {
			__m256i fp_data = _mm256_loadu_si256((__m256i*)&Bid.FP[i]);
			__m256i fp_cmp = _mm256_cmpeq_epi32(fp_data, target_FP);
        	int fp_mask = _mm256_movemask_epi8(fp_cmp);
			if (fp_mask != 0) {
            	match_index = (_tzcnt_u32(fp_mask) / 4)+i;
            	goto search_over;
        	}
			__m256i empty_cmp = _mm256_cmpeq_epi32(fp_data, target_empty);
       		int empty_mask = _mm256_movemask_epi8(empty_cmp);
        	if (empty_mask != 0) {
	 	        first_empty = (_tzcnt_u32(empty_mask) / 4)+i;
				goto search_over; 
		    }
    	}

		for(int index=0;index<S_entry_num;index++){
			if((Bid.p[index]>>2)<p_min_value && (Bid.p[index]&0x1)==0){
				p_min_value = (Bid.p[index]>>2);
				p_min_pos = index;
			}
		}
			

	search_over:
	
		if(match_index>=0){  
			int i = match_index;
			Bid.f[i]++;
			if(Bid.f[i]==pow(2,S_F_len)-1){  
				Bid.f[i]=0;  
				if((Bid.p[i]&0x1)==0){
					Bid.FP[i]=0;
					Bid.f[i]=0;
					Bid.p[i]=0;
					return -1;
				}
				int res = L_vir->report(key,0,0,'f');
				if(res<0){
					Bid.FP[i]=0;
					Bid.f[i]=0;
					Bid.p[i]=0;
					return -1; 
				}				
			}
			if(((Bid.p[i]>>1)&0x1)==0){ 
				uint8_t p_tmp = (Bid.p[i]>>2)+1;
				Bid.p[i] = (Bid.p[i] & 0x03)|((p_tmp & 0x3F)<<2); 
				Bid.p[i] |= 0x2; 
				if((Bid.p[i]>>2)>=P_thr){ 
					if((Bid.p[i]&0x1)==0){ 
						Bid.p[i] |= 0x1; 
						int res = L_vir->report(key,Bid.f[i],Bid.p[i],'o');   
						if(res<0){
							Bid.FP[i]=0;
							Bid.f[i]=0;
							Bid.p[i]=0;
							return -1; 
						}
					}
					else{
						int res = L_vir->report(key,0,0,'p');  
						if(res<0){
							Bid.FP[i]=0;
							Bid.f[i]=0;
							Bid.p[i]=0;
							return -1;  
						}
					}
					Bid.p[i] &= 0x3;  
				}
			}
			return 1; 
		}		
		
	    
		if(first_empty>=0){		
			Bid.FP[first_empty]=keyfp;
			Bid.f[first_empty]=1;
			Bid.p[first_empty]=0;
			Bid.p[first_empty]|=0x06;  
	        return 1;  
		}
		else{ 
			bool kick_op = p_min_value<P_thr;
			static uint32_t rng_state = static_cast<uint32_t>(time(0));
			rng_state ^= rng_state << 13;
        	rng_state ^= rng_state >> 17;
        	rng_state ^= rng_state << 5;
       		uint32_t random_number = rng_state;
			double random_value = random_number / static_cast<double>(UINT32_MAX);
		
			if(LPS_POS){ 
				if (random_value > 1.0/p_min_value) {kick_op=0;}  
			}
			if(kick_op){  
				Bid.FP[p_min_pos] = keyfp;
				Bid.f[p_min_pos]=1;
				Bid.p[p_min_pos]=0;
				Bid.p[p_min_pos]|=0x06;
				return 1;
			}
			return 0; 
		}

	
		return 0;
   }

	void nextlevel(L_virtual * nl){
        L_vir = nl;
    }
	
};





class L_Sketch: public L_virtual{ 
public:
    Bucket_L *PIRecord;
	int cur_ptr;
	int cur_max;
	uint32_t tmp_w;

    L_Sketch(){
        PIRecord = new Bucket_L();		
		cur_ptr=0;
		tmp_w = 0;
    }

#define Burst_Multi 1.1  

int hash_col = 0;
int no_col =0;

    int report(uint32_t ID, uint8_t f,uint8_t p,char mode) {

		int op = -1;
		no_col ++;
	
		if(mode=='o'){ 
			PIRecord->ID[cur_ptr]=ID;
			PIRecord->of[cur_ptr]=0x1;   
			cur_max++;
			while(PIRecord->ID[cur_ptr]!=0)cur_ptr++;
			return 1;  
		}
		else if(mode=='f'){
			for(int i=0;i<cur_max;i++){
				if(PIRecord->ID[i]==ID){    
					if((PIRecord->of[i]>>L_f_of_len)>=pow(2,L_f_of_len)){
						hash_col++;
						goto post;
					}
					uint8_t f_tmp = (PIRecord->of[i] >> L_p_of_len)+1;
					PIRecord->of[i] = (PIRecord->of[i] & 0xFF | ((f_tmp&0xFF)<<L_p_of_len));
					op = i;
					goto post;
				}
			}
			hash_col++;
		}
		else if(mode=='p'){ 
			for(int i=0;i<cur_max;i++){
				if(PIRecord->ID[i]==ID){
					if((PIRecord->of[i]&0xFF)>=pow(2,L_p_of_len)){
						hash_col++;
						goto post;
					}
					PIRecord->of[i] = (PIRecord->of[i] & 0xFF00) | ((PIRecord->of[i] + 1) & 0x00FF);
					op = i;
					goto post;
				}
			}
			hash_col++;
		}
		else{
			assert(0);
		}

	post:
		if(op<0){return 0;}  
		int cur_f = PIRecord->of[op]>>L_f_of_len;
		int cur_p = PIRecord->of[op]&0xFF;
		if(cur_f>=cur_p && cur_p<pow(2,L_p_of_len)){
			PIRecord->ID[op]=0;
			PIRecord->of[op]=0;
			if(cur_ptr>op) cur_ptr = op;   
			return -1;    
		}		

		return 1; 
		
    }

};


class LPSSketch{
public:
	S_Sketch S_layer;
	L_Sketch L_layer;
	BOBHash32 *P_Hash;

	LPSSketch(const S_Sketch& SL, const L_Sketch& LL) : S_layer(SL), L_layer(LL){
		S_layer.nextlevel(&L_layer);
		P_Hash = new BOBHash32(97);
	}

	bool insert(uint32_t key){	
		S_layer.insert(key);
		return 0;
	}		

	vector<report_item> getPI(int mode){
		vector<report_item> res;
		report_item new_PI;
		for(int i=0;i<L_record_num;i++){
			L_ID_type ID = L_layer.PIRecord->ID[i];
			if(ID){
				int f = ((L_layer.PIRecord->of[i])>>L_p_of_len)*(pow(2,S_F_len)-1);
				int p = ((L_layer.PIRecord->of[i])&0xFF)*P_thr;
				S_ID_type FP = ID%pow(2,L1_ID_len);
				int hash_index = P_Hash->run((const char *)&ID, 4) % S_bucket_num;
					
				Bucket_S Bid = S_layer.Filter[hash_index];
				for(int j=0;j<S_entry_num;j++){
					if(Bid.FP[j]==ID){
						f+= Bid.f[j];
						p+= (Bid.p[j]>>2); 
					}
				}
				double den = (double)f/p;
				if(mode){
					if(den<D_thr){ 
						new_PI = {ID,p,den};
						res.push_back(new_PI);
					}
					else{
						;
					}
				}
				else{
					new_PI = {ID,p,den};
					res.push_back(new_PI);
				}
			}
		}
		return res;
	}
	
	bool clear();
	bool init();

};


#endif
