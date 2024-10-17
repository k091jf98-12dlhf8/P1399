  #ifndef CLASS_H
#define CLASS_H

#include "hash.h"
#include "BOBHash32.h"
#include <cstring>
#include <algorithm>
#include <vector>
#include <stdint.h>
#include "para.h"
#include <time.h>
#include <stdlib.h>

using namespace std;

int kick_p = 0;
int kick_t = 0;
int kick_max=-1;

uint32_t getFP(uint32_t key, int key_len)
{
    static BOBHash32 fpHash(100);
    return fpHash.run((const char *)&key, 4) % 0xFFFF + 1;
}



uint32_t cur_w=0; 


struct PSS_RES{
	int flag;   
	L1_F_type f;
	L1_P_type p;
	int index;
	int pos;
};

class Bucket_PSS
{
public:
    L1_ID_type *ID; 
    L1_F_type *f;  
	L1_P_type *p; 
	//uint16_t *w; 
	bool *wflag;
	double *den; 
	uint32_t *kick_time; 
	bool *if_L2;   

    Bucket_PSS() {}

    Bucket_PSS(int entry_num)     
    {
        ID = new L1_ID_type[entry_num];
        f = new L1_F_type[entry_num];
		p = new L1_P_type[entry_num];
		//w = new uint16_t[entry_num];
		wflag = new bool[entry_num];
		den = new double[entry_num];
		kick_time = new uint32_t[entry_num];
		
		if_L2 = new bool[entry_num]();
        memset(ID, 0, sizeof(L1_ID_type) * entry_num);
		memset(f, 0, sizeof(L1_F_type) * entry_num);
		memset(p, 0, sizeof(L1_P_type) * entry_num);
		//memset(w, 0, sizeof(uint16_t) * entry_num);
		memset(wflag, 0, sizeof(bool) * entry_num);
		memset(den, 0, sizeof(double) * entry_num);
		memset(kick_time, 0, sizeof(uint32_t) * entry_num);
    }

};


int lock_b = 0;
int lock[L1_bucket_num]={0};


class L2_virtual_PSS
{
public:
	virtual bool insert(uint32_t ID, int f,int p, bool if_L2) = 0; 
};

class P_Filter{
public:
    Bucket_PSS *Filter;  
    int *kick_time;  
	
    BOBHash32 *P_Hash;
    int bucket_num;
    int entry_num; 
    int P_threshold;  
	int w_cnt;

	L2_virtual_PSS* L2_vir;

    P_Filter() {}
    P_Filter(int _bucket_num, int _entry_num,
                 int _ID_len, int _f_len, int _p_len, int _w_len,
                 int _P_threshold,
                 int rand_seed)
        : bucket_num(_bucket_num), entry_num(_entry_num),
          ID_len(_ID_len), f_len(_f_len),p_len(_p_len),w_len(_w_len),
          P_threshold(_P_threshold)
    {
        Filter = new Bucket_PSS[bucket_num];
		kick_time = new int[bucket_num];
        for (int i = 0; i < bucket_num; ++i){
            Filter[i] = Bucket_PSS(entry_num);
			kick_time[i]=0;
        }
        
        P_Hash = new BOBHash32(rand_seed);
		w_cnt = 0;
    }

	inline bool report(int ID,int f,int p, bool if_L2){ 
		if(if_L2){
			return L2_vir->insert(ID,f,p,if_L2);  
		}  
		return L2_vir->insert(ID,f,p,if_L2);		
	}

	void check_full(int h,int cnt){
		Bucket_PSS Bid = Filter[h];
		bool flag=1;
		for(int i=0;i<entry_num;i++){
			if(Bid.p[i]<P_thr){flag=0;}
		}
		if(flag && lock[h]==0){
			lock_b++;
			lock[h]=1;
		}
	}


	PSS_RES insert(uint32_t key, int count = 1) 
    {
    	w_cnt++;
		if(w_cnt%WindowSize==0){
			cur_w++;
			for(int i=0;i<bucket_num;i++){
				for(int j=0;j<entry_num;j++){
					Filter[i].wflag[j]=0;
				}	
			}
		}
		L1_ID_type keyfp = key%pow(2,L1_ID_len);
		
		PSS_RES re = {-1,0,0,0,0};

        int hash_index[5] = {
        	P_Hash->run((const char *)&key, 4) % bucket_num,
			(P_Hash->run((const char *)&key, 4)+1) % bucket_num,
			(P_Hash->run((const char *)&key, 4)+2) % bucket_num,
			(P_Hash->run((const char *)&key, 4)+3) % bucket_num,
			(P_Hash->run((const char *)&key, 4)+4) % bucket_num
		};
		Bucket_PSS Bid[5] = {
			Filter[hash_index[0]],
			Filter[hash_index[1]],
			Filter[hash_index[2]],
			Filter[hash_index[3]],
			Filter[hash_index[4]]
		};

		for(int h=0;h<5;h++){	

	        for (int i = 0; i < entry_num; ++i)  
	            if (Bid[h].ID[i] == keyfp)  
	            {
	            	Bid[h].f[i]++;
					if(Bid[h].wflag[i]==0){
						Bid[h].p[i]+=1;
						Bid[h].wflag[i]=1;
						if(Bid[h].p[i]>P_threshold){  
							bool report_res = report(key,Bid[h].f[i],Bid[h].p[i],Bid[h].if_L2[i]);  
							Bid[h].if_L2[i] = report_res?1:Bid[h].if_L2[i];  
						}
					}
					re = {1,Bid[h].f[i],Bid[h].p[i],keyfp,i};
					Bid[h].den[i]=(double)(Bid[h].f[i])/Bid[h].p[i];
					check_full(hash_index[h],w_cnt);
					return re;  
	            }
		}

		for(int h=0;h<5;h++){		
	        for (int i = 0; i < entry_num; ++i)
	            if (Bid[h].f[i] == 0)  
	            {
	            	Bid[h].ID[i]=keyfp;
					Bid[h].f[i]=1;
					Bid[h].p[i]=1;
					//Bid.w[i]=cur_w;
					Bid[h].wflag[i]=1;
					
					re = {0,1,1,keyfp,i};
					Bid[h].den[i]=(double)(Bid[h].f[i])/Bid[h].p[i];
	                return re;  
	            }
		}

		uint32_t p_min = (1u<<31)-1;
		int kick[2] = {-1,-1};
		for(int h=0;h<5;h++){		
			for(int i=0;i<entry_num;++i){  
				if(Bid[h].p[i]<p_min && Bid[h].p[i]<P_threshold){
					p_min = Bid[h].p[i];
					kick[0] = h;
					kick[1] = i;
				}
			}
		}
		srand(time(0)); 
		if(kick[0]>=0 && kick[1]>=0){ 
			re = {2,Bid[kick[0]].f[kick[1]],Bid[kick[0]].p[kick[1]],keyfp,kick[1]};
			kick_p += Bid[kick[0]].p[kick[1]];
			kick_t ++;
			if(Bid[kick[0]].p[kick[1]]>kick_max){kick_max=Bid[kick[0]].p[kick[1]];}
			Bid[kick[0]].ID[kick[1]] = keyfp; 
			Bid[kick[0]].f[kick[1]]=1;
			Bid[kick[0]].p[kick[1]]=1; 
			Bid[kick[0]].wflag[kick[1]]=1;
			Bid[kick[0]].kick_time[kick[1]]++; 

			Bid[kick[0]].den[kick[1]]=(double)(Bid[kick[0]].f[kick[1]])/Bid[kick[0]].p[kick[1]];
			kick_time[kick[0]]++;
	        return re;
		}
		
		return re;
		
    }


	void nextlevel(L2_virtual_PSS * nl)
    {
        L2_vir = nl;
    }


	
};



struct Bucket_L2{
    uint32_t fp=0;
    L2_D_type den=0.0; 
    int kick_time=0;
};


class Freq_Sketch: public L2_virtual_PSS{ 
public:
    Bucket_L2 *PS_record;
	int record_num;
	
	int cur_ptr;
	L2_D_type shared_den;

	uint32_t tmp_w;

    Freq_Sketch() {}
    Freq_Sketch(int _record_num)
        : record_num(_record_num)
    {
        PS_record = new Bucket_L2[record_num];		
		cur_ptr=0;
		shared_den=0.0;
		tmp_w = 0;
    }

	vector<uint32_t> getPI_report(){    
		vector<uint32_t> PIList;
		int sum=0;
		for(int i=0;i<cur_ptr;i++){
			if(PS_record[i].den<D_thr){
				PIList.push_back(PS_record[i].fp);
				sum++;
			}
		}
		printf("PI Report:%d  Totally L2:%d\n",sum,cur_ptr);
		return PIList;
	}

	int Rec(){  
		int sum=0;
		for(int i=0;i<cur_ptr;i++){
			if(PS_record[i].den<D_thr){sum++;}
		}
		return sum;
	}

	int ARE(){ 
		return cur_ptr;
	}

	uint32_t getfp(int pos){
		return PS_record[pos].fp;
	}

	double getden(int pos){
		return PS_record[pos].den;
	}

	int getPI(){
		uint32_t res;
		for(int i=0;i<cur_ptr;i++){
			res = PS_record[i].fp;
		}
		return cur_ptr;
	}


#define Burst_Multi 1.1 


    bool insert(uint32_t ID, int f,int p,bool if_kicked) {

		double new_den = (if_kicked&&shared_den>0)?shared_den:((double)f/p);
			
		for(int pos=0;pos<cur_ptr;pos++){ 
			if(PS_record[pos].fp==ID){
				if((double)f/p>Burst_Multi*PS_record[pos].den){
					return 1;
				}
				PS_record[pos].den=(double)f/p;  
				if(PS_record[pos].den==0){printf("?\n");}
				return 1;
			}
		}
	
		int new_pos = cur_ptr; 
		if(cur_ptr==record_num){  
			double den_max = 0.0;
			for(int i=0;i<cur_ptr;i++){
				if(PS_record[i].den>den_max){
					den_max=PS_record[i].den;  
					new_pos = i;
				}
			}
			
			if(new_den >= den_max){return 0;} 
			shared_den += den_max;  
			shared_den /= 2;
			cur_ptr--;
			PS_record[new_pos].kick_time++;
		}
		PS_record[new_pos].fp=ID;
		PS_record[new_pos].den=new_den;

		cur_ptr++;
		if(PS_record[new_pos].den==0){printf("!\n");}
		return 1;
		
    }

};


template<typename L1_Class, typename L2_Class>
class PSSketch{
public:
	P_Filter L1;
	Freq_Sketch L2;

	PSSketch() {}
	PSSketch(const P_Filter& L1_arg, const Freq_Sketch& L2_arg) : L1(L1_arg), L2(L2_arg){
		L1.nextlevel(&L2);
	}

	bool insert(uint32_t key){
		L1.insert(key);
		return 0;
	}
	vector<report_item> getL2PI_report(){   
		vector<report_item> res;
		vector<uint32_t> reportID = L2.getPI_report();
		for(uint32_t item : reportID){
			report_item rep;
			rep.ID = item;
			for(int i=0;i<L1_bucket_num;i++){
				for(int j=0;j<L1_entry_num;j++){
					if(item%pow(2,L1_ID_len)==L1.Filter[i].ID[j]){rep.p=L1.Filter[i].p[j];goto findD;}
				}
			}
		findD:
			for(int k=0;k<L2_recourd_num;k++){
				if(item == L2.PS_record[k].fp){rep.D=L2.PS_record[k].den;break;}
			}
			res.push_back(rep);
		}
		return res;
	}
	
	int query(int op){
		if(op==1){ 
			printf("No use\n");
		}
		else if(op==2){ 
			return L2.Rec();
		}
		else if(op==3){
			return L2.ARE();
		}
		else{
			printf("No such op\n");
		}
		return 0;
	}
	int getL1f(uint32_t key){
		for(int i=0;i<L1.bucket_num;i++){
			for(int j=0;j<L1.entry_num;j++){
				if(key%pow(2,L1_ID_len)==L1.Filter[i].ID[j]){return L1.Filter[i].f[j];}
			}
		}
		return -1;
	}
	int getL1p(uint32_t key){
		for(int i=0;i<L1.bucket_num;i++){
			for(int j=0;j<L1.entry_num;j++){
				if(key%pow(2,L1_ID_len)==L1.Filter[i].ID[j]){return L1.Filter[i].p[j];}
			}
		}
		return -1;
	}
	uint32_t getL2fp(int pos){
		return L2.getfp(pos);
	}
	double getL2den(int pos){
		return L2.getden(pos);
	}
	int getPI(){
		return L2.getPI();
	}
	
	bool clear();
	bool init();

};






#endif
