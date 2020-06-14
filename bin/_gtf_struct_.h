#ifndef _GTF_STRUCT_H
#define _GTF_STRUCT_H

#include <string>
#include <map>
#include <iostream>
#include <memory>

#include <unordered_map>
#include <unordered_set>
#include <set>

#include "_com_structs_.h"

//enum _gtf_entry_type_ {"gene","transcript","exon","CDS","UTR","start_codon","stop_codon","Selenocysteine"};

class _gtf_entry_{
protected:
   std::string _this_type;
   int64_t _coding_length;

   std::string _name_;
   std::string _type_;
   std::string _id_;

   GenomicRegion _this_gr;

   //std::unordered_map<int64_t, bool>::iterator _check_pos_it;
   //std::unordered_set<int64_t>::iterator _check_pos_it;
   std::set<int64_t>::iterator _check_pos_it;

public:
   std::vector< std::shared_ptr<_gtf_entry_> > _sub_list;

   int64_t get_length(){ return _this_gr.start_pos - _this_gr.end_pos; }
   std::string get_type() {return _this_type; }
   void set_name(std::string _p_name){ _name_ = _p_name; }
   void set_gene_type(std::string _p_type){ _type_ = _p_type; }
   void set_id(std::string _p_id){ _id_ = _p_id; }
   std::string get_name(){ return _name_; }
   std::string get_gene_type(){ return _type_; }
   std::string get_id(){ return _id_; }
   void set_this_gr(const GenomicRegion _p_this_gr){ _this_gr = _p_this_gr; }
   GenomicRegion get_this_gr() { return _this_gr; }

   void set_gr_start_pos(const int64_t _p_pos);
   void set_gr_end_pos(const int64_t _p_pos);
   void set_gr_chr(const std::string _p_chrn){ _this_gr.chrn = _p_chrn; }
   int64_t get_gr_start_pos(){ return  _this_gr.start_pos; }
   int64_t get_gr_end_pos(){ return  _this_gr.end_pos; }
   std::string get_gr_chr(){ return  _this_gr.chrn; }

   virtual int64_t get_coding_length(){ return 100;}
   virtual int64_t get_coding_ovlp(GenomicMapRegion & _gmr_){return 100;}
   //virtual int64_t get_coding_ovlp(std::map<int64_t, bool>& _this_map_info, std::string _t_chrn){return 100;}
   //virtual int64_t get_coding_ovlp(std::unordered_map<int64_t, bool>& _this_map_info, std::string _t_chrn){return 100;}  
   //virtual int64_t get_coding_ovlp(std::unordered_set<int64_t>& _this_map_info, std::string _t_chrn){return 100;}
   virtual int64_t get_coding_ovlp(std::set<int64_t>& _this_map_info, std::string _t_chrn){return 100;}

   _gtf_entry_() { _coding_length = -1; _this_gr.start_pos=0; _this_gr.end_pos=0; } 
   _gtf_entry_(std::string _t_type, std::string _p_n, std::string _p_t, std::string _p_id, GenomicRegion _p_gr);
   _gtf_entry_(std::string _t_type, std::string _p_id) { _this_type = _t_type; set_id(_p_id); _coding_length = -1; _this_gr.start_pos=0; _this_gr.end_pos=0; }
   int64_t get_ovlp(GenomicMapRegion & _gmr_);
   int64_t get_ovlp(std::shared_ptr<_gtf_entry_>  _other_e);
   std::string _to_string();
};

class _gene_entry_ : public _gtf_entry_ {
public:
   int64_t get_coding_length();
   int64_t get_coding_ovlp(GenomicMapRegion & _gmr_);
   //int64_t get_coding_ovlp(std::map<int64_t, bool>& _this_map_info, std::string _t_chrn);
   //int64_t get_coding_ovlp(std::unordered_map<int64_t, bool>& _this_map_info, std::string _t_chrn);
   //int64_t get_coding_ovlp(std::unordered_set<int64_t>& _this_map_info, std::string _t_chrn);
   int64_t get_coding_ovlp(std::set<int64_t>& _this_map_info, std::string _t_chrn);

   _gene_entry_(std::string _t_type, std::string _p_n, std::string _p_t, std::string _p_id, GenomicRegion _p_gr):_gtf_entry_(_t_type,_p_n,_p_t,_p_id,_p_gr){;}
   _gene_entry_(std::string _t_type, std::string _p_id):_gtf_entry_(_t_type,_p_id){ _coding_length = -1;}
};


class _transcript_entry_ : public _gtf_entry_ {
public:
   int64_t get_coding_length();
   int64_t get_coding_ovlp(GenomicMapRegion & _gmr_);
   //int64_t get_coding_ovlp(std::map<int64_t, bool>& _this_map_info, std::string _t_chrn);
   //int64_t get_coding_ovlp(std::unordered_map<int64_t, bool>& _this_map_info, std::string _t_chrn);
   //int64_t get_coding_ovlp(std::unordered_set<int64_t>& _this_map_info, std::string _t_chrn);
   int64_t get_coding_ovlp(std::set<int64_t>& _this_map_info, std::string _t_chrn);

   _transcript_entry_(std::string _t_type, std::string _p_n, std::string _p_t, std::string _p_id, GenomicRegion _p_gr):_gtf_entry_(_t_type,_p_n,_p_t,_p_id,_p_gr){;}
   _transcript_entry_(std::string _t_type, std::string _p_id):_gtf_entry_(_t_type,_p_id){ _coding_length = -1;}
};

class _exon_entry_ : public _gtf_entry_ {
public:
   int64_t get_coding_length();
   int64_t get_coding_ovlp(GenomicMapRegion & _gmr_);
   //int64_t get_coding_ovlp(std::map<int64_t, bool>& _this_map_info, std::string _t_chrn);
   //int64_t get_coding_ovlp(std::unordered_map<int64_t, bool>& _this_map_info, std::string _t_chrn);
   //int64_t get_coding_ovlp(std::unordered_set<int64_t>& _this_map_info, std::string _t_chrn);
   int64_t get_coding_ovlp(std::set<int64_t>& _this_map_info, std::string _t_chrn);

   _exon_entry_(std::string _t_type, std::string _p_n, std::string _p_t, std::string _p_id, GenomicRegion _p_gr):_gtf_entry_(_t_type,_p_n,_p_t,_p_id,_p_gr){;}
   //_exon_entry_(std::string _t_type):_gtf_entry_(_t_type){ _coding_length = -1;}
}; 

class _other_entry_ : public _exon_entry_ {
public:
   int64_t get_coding_length();
   int64_t get_coding_ovlp(GenomicMapRegion & _gmr_);
   //int64_t get_coding_ovlp(std::map<int64_t, bool>& _this_map_info, std::string _t_chrn);
   //int64_t get_coding_ovlp(std::unordered_map<int64_t, bool>& _this_map_info, std::string _t_chrn);
   //int64_t get_coding_ovlp(std::unordered_set<int64_t>& _this_map_info, std::string _t_chrn);
   int64_t get_coding_ovlp(std::set<int64_t>& _this_map_info, std::string _t_chrn);

   _other_entry_(std::string _t_type, std::string _p_n, std::string _p_t, std::string _p_id, GenomicRegion _p_gr):_exon_entry_(_t_type,_p_n,_p_t,_p_id,_p_gr){;}
   //_other_entry_(std::string _t_type):_exon_entry_(_t_type){ _coding_length = -1;}
};


int read_gtf(std::map<std::string, std::map<std::string,  std::shared_ptr<_gtf_entry_> > > & _gene_dict, std::map<std::string, std::string> & _gid_to_gn, std::string _in_gtf_file, const char ** gt_list, const int gt_size);


#endif

