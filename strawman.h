
#include "CMSketch.h"
#include "OO_PE.h"
#include <vector>
#include "para.h"

class Strawman {
public:
	CMSketch cmSketch;
    OO_PE<uint32_t, int16_t> oo_pe;
    vector<pair<uint32_t,double>> PI_List;
	int w_cnt;
	int w_num;

    Strawman()
    : cmSketch(OO_PE_Len, 3), oo_pe(3, OO_PE_Len) {
		w_cnt=0;
		w_num=0;
	}

    void insert(uint32_t item) {
    	w_cnt++;
		if(w_cnt%WindowSize==0){
			w_num++;
			oo_pe.NewWindow(w_num); 
		}
		int p = oo_pe.Insert(item,w_num);
        pair<uint32_t,double> res = cmSketch.Insert(item,p);
		if(res.second>0){ 
			bool found = false;
        	for (auto& find : PI_List) {
            if (find.first == res.first) { 
                find.second = res.second;
                found = true;
                break;
            	}
        	}
        	if (!found) {
            	PI_List.push_back({res.first, res.second});
        	}
		} 
    }

   uint32_t query(uint32_t item) {
        uint32_t oo_pe_count = oo_pe.Query(item);
        int cm_count = cmSketch.Query(reinterpret_cast<const char*>(&item));
        return min(static_cast<uint32_t>(cm_count), oo_pe_count);
   }

   void newWindow() {
        oo_pe.NewWindow(0);
   }

    
};




































