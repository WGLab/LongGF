
#include "_gtf_struct_.h"
#include "_com_fun_.h"

#include <sstream>
#include <fstream>
#include <string>
#include <iostream>


_gtf_entry_::_gtf_entry_(std::string _t_type, std::string _p_n, std::string _p_t, std::string _p_id, GenomicRegion _p_gr){
      _this_type = _t_type;
      _coding_length = -1;
      set_name(_p_n);
      set_gene_type(_p_t);
      set_id(_p_id);
      set_this_gr(_p_gr);
}

void _gtf_entry_::set_gr_start_pos(const int64_t _p_pos){
   if ((_this_gr.start_pos==0 || _p_pos<_this_gr.start_pos) && _p_pos>=0){
      _this_gr.start_pos = _p_pos;
   }
}
void _gtf_entry_::set_gr_end_pos(const int64_t _p_pos){
   if ((_this_gr.end_pos==0 || _p_pos>_this_gr.end_pos) && _p_pos>=0){
      _this_gr.end_pos = _p_pos;
   }
}

std::string _gtf_entry_::_to_string(){
   return _this_type+" "+_name_+" "+_type_+" "+_id_+" "+std::to_string(_coding_length)+"/"+std::to_string(_sub_list.size())+" "+_this_gr.chrn+":"+std::to_string(_this_gr.start_pos)+"-"+std::to_string(_this_gr.end_pos);
}

int64_t _gtf_entry_::get_ovlp(GenomicMapRegion & _gmr_){
   //std::cout<<_gmr_.chrn<<":"<<_gmr_.ref_start_pos<<":"<<_gmr_.ref_end_pos<< " "<<_this_gr.chrn<<":"+std::to_string(_this_gr.start_pos)+"-"+std::to_string(_this_gr.end_pos)<<" ";
   if (_gmr_.chrn.compare(_this_gr.chrn)==0){
      int64_t max_start = (_gmr_.ref_start_pos>_this_gr.start_pos?_gmr_.ref_start_pos:_this_gr.start_pos);
      int64_t min_end   = (_gmr_.ref_end_pos  <_this_gr.end_pos  ?_gmr_.ref_end_pos  :_this_gr.end_pos);
      //std::cout<<min_end-max_start<<std::endl;
      if (min_end-max_start>0){ return min_end-max_start; }
      else{ return 0; }
   }else{
      //std::cout<<0<<std::endl;
      return 0;
   }
}
int64_t _gtf_entry_::get_ovlp(std::shared_ptr<_gtf_entry_>  _other_e){
   if (_other_e->_this_gr.chrn.compare(_this_gr.chrn)==0){
      int64_t max_start = (_other_e->_this_gr.start_pos>_this_gr.start_pos?_other_e->_this_gr.start_pos:_this_gr.start_pos);
      int64_t min_end   = (_other_e->_this_gr.end_pos  <_this_gr.end_pos  ?_other_e->_this_gr.end_pos  :_this_gr.end_pos);
      if (min_end-max_start>0){ return min_end-max_start; }
      else{ return 0; }
   }else{
      return 0;
   }
}

int64_t _gene_entry_::get_coding_length(){
   if (_coding_length>0) { return _coding_length; }
   else{
      _coding_length = 0;
      int64_t _t_len = 0;
      for (std::vector< std::shared_ptr<_gtf_entry_> >::iterator _f_it = _sub_list.begin(); _f_it != _sub_list.end(); _f_it++){
         _t_len = (*_f_it)->get_coding_length();
         if (_t_len>_coding_length){
             _coding_length = _t_len;
         }
      }
      return _coding_length;
   }
}
int64_t _gene_entry_::get_coding_ovlp(GenomicMapRegion & _gmr_){
   //std::cout<<"Testovlp g="<<_gmr_.chrn<<":"<<_gmr_.ref_start_pos<<":"<<_gmr_.ref_end_pos<< " "<<_this_gr.chrn<<":"+std::to_string(_this_gr.start_pos)+"-"+std::to_string(_this_gr.end_pos)<<" "<<_gmr_.chrn.compare(_this_gr.chrn)<<std::endl;
   if (_gmr_.chrn.compare(_this_gr.chrn)==0){
      int64_t _t_len = 0;
      for (std::vector< std::shared_ptr<_gtf_entry_> >::iterator _f_it = _sub_list.begin(); _f_it != _sub_list.end(); _f_it++){
         int64_t _c_len =(*_f_it)->get_coding_ovlp(_gmr_);
         if (_c_len>_t_len){
             _t_len = _c_len;
         }
      }
      return _t_len;
   }else{
      return 0;
   }
}
//int64_t _gene_entry_::get_coding_ovlp(std::map<int64_t, bool>& _this_map_info, std::string _t_chrn){
//int64_t _gene_entry_::get_coding_ovlp(std::unordered_map<int64_t, bool>& _this_map_info, std::string _t_chrn){
//int64_t _gene_entry_::get_coding_ovlp(std::unordered_set<int64_t>& _this_map_info, std::string _t_chrn){
int64_t _gene_entry_::get_coding_ovlp(std::set<int64_t>& _this_map_info, std::string _t_chrn){
  if (_t_chrn.compare(_this_gr.chrn)==0){
      int64_t _t_len = 0;
      for (std::vector< std::shared_ptr<_gtf_entry_> >::iterator _f_it = _sub_list.begin(); _f_it != _sub_list.end(); _f_it++){
         int64_t _c_len =(*_f_it)->get_coding_ovlp(_this_map_info, _t_chrn);
         if (_c_len>_t_len){
             _t_len = _c_len;
         }
      }
      return _t_len;
   }else{
      return 0;
   } 
}


int64_t _exon_entry_::get_coding_length(){
   int64_t _e_l = get_length();
   if (_e_l>0){ return _e_l; }
   else{ return 0; }
}
int64_t _exon_entry_::get_coding_ovlp(GenomicMapRegion & _gmr_){
   //std::cout<<"\t\tTestovlp e=";
   return get_ovlp(_gmr_); 

   //int64_t max_start = (_gmr_.ref_start_pos>_this_gr.start_pos?_gmr_.ref_start_pos:_this_gr.start_pos);
   //int64_t min_end   = (_gmr_.ref_end_pos  <_this_gr.end_pos  ?_gmr_.ref_end_pos  :_this_gr.end_pos);
   //if (min_end-max_start>0){ return min_end-max_start; }
   //else{ return 0; }
}
//int64_t _exon_entry_::get_coding_ovlp(std::map<int64_t, bool>& _this_map_info, std::string _t_chrn){
//int64_t _exon_entry_::get_coding_ovlp(std::unordered_map<int64_t, bool>& _this_map_info, std::string _t_chrn){
//int64_t _exon_entry_::get_coding_ovlp(std::unordered_set<int64_t>& _this_map_info, std::string _t_chrn){
int64_t _exon_entry_::get_coding_ovlp(std::set<int64_t>& _this_map_info, std::string _t_chrn){
   int64_t _t_ovlp_l = 0;
  
   for (int64_t _t_pos=_this_gr.start_pos-1; _t_pos<_this_gr.end_pos; _t_pos++){
      _check_pos_it = _this_map_info.find(_t_pos);
      if (_check_pos_it != _this_map_info.end()){
          _t_ovlp_l += 1;
      }
   }
 
   return _t_ovlp_l; 
}

int64_t _transcript_entry_::get_coding_length(){
   if (_coding_length>0) { return _coding_length; }
   else{
      _coding_length = 0;
      for (std::vector< std::shared_ptr<_gtf_entry_> >::iterator _f_it = _sub_list.begin(); _f_it != _sub_list.end(); _f_it++){
         _coding_length += (*_f_it)->get_coding_length();
      }
      return _coding_length;
   }
}
int64_t _transcript_entry_::get_coding_ovlp(GenomicMapRegion & _gmr_){
   int64_t _t_len = 0;
   for (std::vector< std::shared_ptr<_gtf_entry_> >::iterator _f_it = _sub_list.begin(); _f_it != _sub_list.end(); _f_it++){
      _t_len += (*_f_it)->get_coding_ovlp(_gmr_);
   }
   //std::cout<<"\tTestovlp t="<<_gmr_.chrn<<":"<<_gmr_.ref_start_pos<<":"<<_gmr_.ref_end_pos<< " "<<_this_gr.chrn<<":"+std::to_string(_this_gr.start_pos)+"-"+std::to_string(_this_gr.end_pos)<<" " <<_t_len<<std::endl;
   return _t_len;
}
//int64_t _transcript_entry_::get_coding_ovlp(std::map<int64_t, bool>& _this_map_info, std::string _t_chrn){
//int64_t _transcript_entry_::get_coding_ovlp(std::unordered_map<int64_t, bool>& _this_map_info, std::string _t_chrn){
//int64_t _transcript_entry_::get_coding_ovlp(std::unordered_set<int64_t>& _this_map_info, std::string _t_chrn){
int64_t _transcript_entry_::get_coding_ovlp(std::set<int64_t>& _this_map_info, std::string _t_chrn){
   int64_t _t_len = 0;
   for (std::vector< std::shared_ptr<_gtf_entry_> >::iterator _f_it = _sub_list.begin(); _f_it != _sub_list.end(); _f_it++){
      _t_len += (*_f_it)->get_coding_ovlp(_this_map_info, _t_chrn);
   }
   return _t_len;
}

int64_t _other_entry_::get_coding_length(){
   return 0;
   
   int64_t _e_l = get_length();
   if (_e_l>0){ return _e_l; }
   else{ return 0; }
}
int64_t _other_entry_::get_coding_ovlp(GenomicMapRegion & _gmr_){
   return 0;
}
//int64_t _other_entry_::get_coding_ovlp(std::map<int64_t, bool>& _this_map_info, std::string _t_chrn){
//int64_t _other_entry_::get_coding_ovlp(std::unordered_map<int64_t, bool>& _this_map_info, std::string _t_chrn){
//int64_t _other_entry_::get_coding_ovlp(std::unordered_set<int64_t>& _this_map_info, std::string _t_chrn){
int64_t _other_entry_::get_coding_ovlp(std::set<int64_t>& _this_map_info, std::string _t_chrn){
   return 0;
}


int read_gtf(std::map<std::string, std::map<std::string,  std::shared_ptr<_gtf_entry_> > > & _gene_dict, std::map<std::string, std::string> & _gid_to_gn, std::string _in_gtf_file, const char ** gt_list, const int gt_size){
   std::map<std::string, std::map<std::string,  std::shared_ptr<_gtf_entry_> > >::iterator g_d_it;

   std::string line;
   std::ifstream infile(_in_gtf_file);
   std::vector<std::string> m_substr_list;

   //const char * gt_list1[100] = {"IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene","TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene","IG_pseudogene","IG_C_pseudogene","IG_J_pseudogene","IG_V_pseudogene","TR_V_pseudogene","TR_J_pseudogene","nonsense_mediated_decay","non_stop_decay","protein_coding","ambiguous_orf","pseudogene","processed_pseudogene","polymorphic_pseudogene","transcribed_processed_pseudogene","transcribed_unprocessed_pseudogene","transcribed_unitary_pseudogene","translated_processed_pseudogene","translated_unprocessed_pseudogene","unitary_pseudogene","unprocessed_pseudogene","disrupted_domain"};
   //int gt_size1 = 30;
   //
   //const char * gt_list2[100] = {"IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene","TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene","nonsense_mediated_decay","non_stop_decay","protein_coding","ambiguous_orf","disrupted_domain"};
   //int gt_size2 = 14;
   //int gt_size = gt_size2;
   //const char ** gt_list;
   //if (_used_pseudogene==1){
   //   gt_size = gt_size1;
   //   gt_list = gt_list1;
   //}else{
   //   gt_size = gt_size2;
   //   gt_list = gt_list2;
   //}

   std::string cur_g_str;
   while (std::getline(infile, line)){
      if (line.size()>0 && line[0]!='#'){
          m_substr_list = m_split_string(line, '\t');
          if (m_substr_list.size()>0){
             //if (m_substr_list[0].size()<3 || m_substr_list[0].substr(0,3).compare("chr")!=0) {continue; }
             // revise the line above to tolerate chrzz with zz.
             if (m_substr_list[0].size()>=3 && str_tolower(m_substr_list[0].substr(0,3)).compare("chr")==0){
                m_substr_list[0] = m_substr_list[0].substr(3);
             }
             if (m_substr_list[0].find_first_of('_')!=std::string::npos){continue;}

             GenomicRegion cur_gr;
             cur_gr.chrn = m_substr_list[0];
             cur_gr.start_pos = std::stol(m_substr_list[3]);
             cur_gr.end_pos = std::stol(m_substr_list[4]);
             if (cur_gr.end_pos < cur_gr.start_pos ){
                 int64_t tmp = cur_gr.end_pos;
                 cur_gr.end_pos = cur_gr.start_pos;
                 cur_gr.start_pos = tmp;
             }

             std::stringstream oss(m_substr_list[8]);
             std::string cur_g_name;
             std::string cur_g_type;
             std::string cur_g_id;
             std::string cur_t_id;
             std::string cur_t_type;
             std::string cur_t_name;
             std::string cur_e_id;
             std::string cur_e_name;
             cur_g_str.clear();
             while (std::getline(oss, cur_g_str, ' ')){
                 cur_g_str = str_tolower(cur_g_str); // revise to tolerate lowcase and uppercase
                 if (cur_g_str.compare("gene_id")==0){
                     std::getline(oss, cur_g_id, ' ');
                     cur_g_id = cur_g_id.substr(1, cur_g_id.size()-3);
                 }else if (cur_g_str.compare("transcript_id")==0){
                     std::getline(oss, cur_t_id, ' ');
                     cur_t_id = cur_t_id.substr(1, cur_t_id.size()-3);
                 }else if (cur_g_str.compare("gene_name")==0){
                     std::getline(oss, cur_g_name, ' ');
                     cur_g_name = cur_g_name.substr(1, cur_g_name.size()-3);
                 }else if (cur_g_str.compare("gene_type")==0 || cur_g_str.compare("gene_biotype")==0 ){// add gene_biotype to consider other types
                     std::getline(oss, cur_g_type, ' ');
                     cur_g_type = cur_g_type.substr(1, cur_g_type.size()-3);
                 }else if (cur_g_str.compare("transcript_name")==0){
                     std::getline(oss, cur_t_name, ' ');
                     cur_t_name = cur_t_name.substr(1, cur_t_name.size()-3);
                 }else if (cur_g_str.compare("transcript_type")==0){
                     std::getline(oss, cur_t_type, ' ');
                     cur_t_type = cur_t_type.substr(1, cur_t_type.size()-3);
                 }else if (cur_g_str.compare("exon_number")==0){
                     std::getline(oss, cur_e_name, ' ');
                     cur_e_name = cur_e_name.substr(0, cur_e_name.size()-1);
                 }else if (cur_g_str.compare("exon_id")==0){
                     std::getline(oss, cur_e_id, ' ');
                     cur_e_id = cur_e_id.substr(1, cur_e_id.size()-3);
                 }
             }

             bool is_g=false;
             for (int isg_i=0; isg_i<gt_size; isg_i++){
                if (cur_g_type.compare(gt_list[isg_i])==0){
                   is_g = true;
                   break;
                }
             }
             if (gt_size>0){ 
                if (!is_g){continue;}
             }

             if (m_substr_list[2].compare("gene")==0){
                std::shared_ptr<_gene_entry_> _c_g_ = std::make_shared<_gene_entry_>("gene", cur_g_name, cur_g_type, cur_g_id, cur_gr);
                _gid_to_gn[cur_g_id] = (cur_g_name.compare("")==0?cur_g_id:cur_g_name);
                g_d_it = _gene_dict.find(cur_gr.chrn);
                if (g_d_it==_gene_dict.end()){
                    _gene_dict[cur_gr.chrn] = std::map<std::string,  std::shared_ptr<_gtf_entry_> >();
                }
                if (_gene_dict[cur_gr.chrn].find(cur_g_id)==_gene_dict[cur_gr.chrn].end()){
                    _gene_dict[cur_gr.chrn][cur_g_id] = _c_g_;
                }else{
                    _gene_dict[cur_gr.chrn][cur_g_id]->set_name(cur_g_name);
                    _gene_dict[cur_gr.chrn][cur_g_id]->set_gene_type(cur_g_type);
                    _gene_dict[cur_gr.chrn][cur_g_id]->set_id(cur_g_id);
                    _gene_dict[cur_gr.chrn][cur_g_id]->set_this_gr(cur_gr);
                }
             } else if (m_substr_list[2].compare("transcript")==0){
                std::shared_ptr<_transcript_entry_> _c_t = std::make_shared<_transcript_entry_>("transcript", cur_t_name, cur_t_type, cur_t_id, cur_gr);
                g_d_it = _gene_dict.find(cur_gr.chrn);
                if (g_d_it==_gene_dict.end()){
                    _gene_dict[cur_gr.chrn] = std::map<std::string,  std::shared_ptr<_gtf_entry_> >();
                }
                if (_gene_dict[cur_gr.chrn].find(cur_g_id)==_gene_dict[cur_gr.chrn].end()){
                    _gene_dict[cur_gr.chrn][cur_g_id] = std::make_shared<_gene_entry_>("gene", cur_g_id);
                    _gene_dict[cur_gr.chrn][cur_g_id]->set_gr_start_pos(cur_gr.start_pos);
                    _gene_dict[cur_gr.chrn][cur_g_id]->set_gr_end_pos(cur_gr.end_pos);
                    _gene_dict[cur_gr.chrn][cur_g_id]->set_gr_chr(cur_gr.chrn);
                    //std::cout<< "gen="<<_gene_dict[cur_gr.chrn][cur_g_id]->_to_string()<<" " <<_c_t->_to_string()<<std::endl;
                }
                std::vector< std::shared_ptr<_gtf_entry_> >::iterator _sl_it=_gene_dict[cur_gr.chrn][cur_g_id]->_sub_list.begin(); 
                for(; _sl_it!=_gene_dict[cur_gr.chrn][cur_g_id]->_sub_list.end(); _sl_it++){
                    if ((*_sl_it)->get_id().compare(cur_t_id)==0){
                        break;
                    }
                }
                if (_sl_it==_gene_dict[cur_gr.chrn][cur_g_id]->_sub_list.end()){
                    _gene_dict[cur_gr.chrn][cur_g_id]->_sub_list.push_back(_c_t);
                }else{
                    (*_sl_it)->set_name(cur_t_name);
                    (*_sl_it)->set_gene_type(cur_t_type);
                    (*_sl_it)->set_id(cur_t_id);
                    (*_sl_it)->set_this_gr(cur_gr);
                }
             } else if (m_substr_list[2].compare("exon")==0){
                 std::shared_ptr<_exon_entry_> _c_e = std::make_shared<_exon_entry_>("exon", cur_e_name, "", cur_e_id, cur_gr);
                g_d_it = _gene_dict.find(cur_gr.chrn);
                if (g_d_it==_gene_dict.end()){
                    _gene_dict[cur_gr.chrn] = std::map<std::string,  std::shared_ptr<_gtf_entry_> >();
                }
                if (_gene_dict[cur_gr.chrn].find(cur_g_id)==_gene_dict[cur_gr.chrn].end()){
                    _gene_dict[cur_gr.chrn][cur_g_id] = std::make_shared<_gene_entry_>("gene", cur_g_id);
                    _gene_dict[cur_gr.chrn][cur_g_id]->set_gr_start_pos(cur_gr.start_pos);
                    _gene_dict[cur_gr.chrn][cur_g_id]->set_gr_end_pos(cur_gr.end_pos);
                    _gene_dict[cur_gr.chrn][cur_g_id]->set_gr_chr(cur_gr.chrn);
                }
                std::vector< std::shared_ptr<_gtf_entry_> >::iterator _sl_it=_gene_dict[cur_gr.chrn][cur_g_id]->_sub_list.begin();
                for(; _sl_it!=_gene_dict[cur_gr.chrn][cur_g_id]->_sub_list.end(); _sl_it++){
                    if ((*_sl_it)->get_id().compare(cur_t_id)==0){
                        break;
                    }
                }
                if (_sl_it==_gene_dict[cur_gr.chrn][cur_g_id]->_sub_list.end()){
                   std::shared_ptr<_transcript_entry_> _n_t=std::make_shared<_transcript_entry_>("transcript", cur_t_id);
                   _n_t->set_gr_start_pos(cur_gr.start_pos);
                   _n_t->set_gr_end_pos(cur_gr.end_pos);
                   _n_t->set_gr_chr(cur_gr.chrn);

                   _n_t->_sub_list.push_back(_c_e);
                   _gene_dict[cur_gr.chrn][cur_g_id]->_sub_list.push_back(_n_t);
                }else{
                   (*_sl_it)->_sub_list.push_back(_c_e);
                }
             } else {// 
                 std::shared_ptr<_other_entry_> _c_o = std::make_shared<_other_entry_>(m_substr_list[2], cur_e_name, "", cur_e_id, cur_gr);
                g_d_it = _gene_dict.find(cur_gr.chrn);
                if (g_d_it==_gene_dict.end()){
                    _gene_dict[cur_gr.chrn] = std::map<std::string,  std::shared_ptr<_gtf_entry_> >();
                }
                if (_gene_dict[cur_gr.chrn].find(cur_g_id)==_gene_dict[cur_gr.chrn].end()){
                    _gene_dict[cur_gr.chrn][cur_g_id] = std::make_shared<_gene_entry_>("gene", cur_g_id);
                    _gene_dict[cur_gr.chrn][cur_g_id]->set_gr_start_pos(cur_gr.start_pos);
                    _gene_dict[cur_gr.chrn][cur_g_id]->set_gr_end_pos(cur_gr.end_pos);
                    _gene_dict[cur_gr.chrn][cur_g_id]->set_gr_chr(cur_gr.chrn);
                }
                std::vector< std::shared_ptr<_gtf_entry_> >::iterator _sl_it=_gene_dict[cur_gr.chrn][cur_g_id]->_sub_list.begin();
                for(; _sl_it!=_gene_dict[cur_gr.chrn][cur_g_id]->_sub_list.end(); _sl_it++){
                    if ((*_sl_it)->get_id().compare(cur_t_id)==0){
                        break;
                    }
                }
                if (_sl_it==_gene_dict[cur_gr.chrn][cur_g_id]->_sub_list.end()){
                   std::shared_ptr<_transcript_entry_> _n_t=std::make_shared<_transcript_entry_>("transcript", cur_t_id);
                   _n_t->set_gr_start_pos(cur_gr.start_pos);
                   _n_t->set_gr_end_pos(cur_gr.end_pos);
                   _n_t->set_gr_chr(cur_gr.chrn);

                   _n_t->_sub_list.push_back(_c_o);
                   _gene_dict[cur_gr.chrn][cur_g_id]->_sub_list.push_back(_n_t);
                }else{
                   (*_sl_it)->_sub_list.push_back(_c_o);
                }
             }
          }
      }
   }
   infile.close();

   return 0;
}


