#ifndef PARA_H
#define PARA_H

#include <assert.h>
#include <time.h>


unordered_map<uint32_t, int> ground_truth_f;
unordered_map<uint32_t, int> ground_truth_p;
unordered_map<uint32_t, int> ground_truth_w;


#define L1_bucket_num 3000   
#define L1_entry_num 3    
#define L2_recourd_num 10000   

#define L1_ID_len 16
#define L1_ID_type uint16_t
#define L1_F_len 16
#define L1_F_type uint16_t
#define L1_P_len 16
#define L1_P_type uint16_t
#define L1_w_len 1
#define L1_flag_len 1
#define L2_D_len 32
#define L2_D_type float



#define S_bucket_num 240   
#define S_entry_num 32   
#define S_total 501*32
#define L_record_num 5000   

#define S_ID_len 32
#define S_ID_type uint32_t
#define S_F_len 8
#define S_F_type uint8_t
#define S_P_len 6
#define S_P_type uint8_t
#define S_flag_len 2  

#define L_ID_len 32
#define L_ID_type uint32_t
#define L_f_of_len 8
#define L_p_of_len 8  
#define L_of_type uint16_t 

#define S_size 6
#define L_size 6

#define LPS_debug 0
#define LPS_POS 0
#define SIMD 1




#define PI_X 1500   
#define PI_Y 5  
#define BF_len 1000   
#define L 10
#define W_thr 200



#define OO_PE_Len 100000


#define P_thr 50
#define D_thr 1.5
#define WindowSize 1000
bool w_bit;


#define Seed 99
#define TopN 1000
#define MAX_INSERT_PACKAGE 32000000
#define filename "./MAWI.dat"  


#define L1_size 6.25
#define L2_size 8
#define BF_size 1
#define PI_size 10.125
#define CM_size 4
#define OO_size 6.125


#define fmax 17057
#define pmax 2452


inline long pow(int a, int x){ 
	return 1L<<x;
}

struct report_item{
	uint32_t ID;
	int p;
	double D;
};



struct timespec pin_s, pin_e;
double Tsum=0;
int Tcnt=0;
inline void PIN(char p){
	if(p=='s')clock_gettime(CLOCK_MONOTONIC, &pin_s);
	else if(p=='e'){
		clock_gettime(CLOCK_MONOTONIC, &pin_e);
		long seconds = pin_e.tv_sec - pin_s.tv_sec;
		long nanoseconds = pin_e.tv_nsec - pin_s.tv_nsec;
		Tsum += (seconds + nanoseconds*1e-9);
		Tcnt ++;
	}
	else{  
		;
	}
}




#endif



