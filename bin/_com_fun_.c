#include "_com_fun_.h"
#include "_com_structs_.h"

#include <sstream>
#include <fstream>
#include <iostream>

#include <algorithm>
#include <cctype>


char m_cigar_str[] = "MIDNSHP=XB";


std::string str_tolower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); } // correct
                  );
    return s;
}


bool compare_support(const NumStrPair & snp1, const NumStrPair & snp2){
   return snp1._num_<snp2._num_;
}

std::vector<std::string> m_split_string(const std::string & m_str, char m_delimiter){
    std::vector<std::string> m_substr_list;

    std::stringstream oss(m_str);
    std::string m_substr;
    while (std::getline(oss, m_substr, m_delimiter)){
        m_substr_list.push_back(m_substr);
    }

    return m_substr_list;
}

std::map<std::string, std::map<std::string, GenomicRegion> > get_gene_from_gtf(const char* in_gtf_file, const int _used_pseudogene){
   std::map<std::string, std::map<std::string, GenomicRegion> > m_gene_list;
   std::map<std::string, std::map<std::string, GenomicRegion> >::iterator gl_it;

   std::string line;
   std::ifstream infile(in_gtf_file);
   std::vector<std::string> m_substr_list;

   const char * gt_list1[100] = {"IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene","TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene","IG_pseudogene","IG_C_pseudogene","IG_J_pseudogene","IG_V_pseudogene","TR_V_pseudogene","TR_J_pseudogene","nonsense_mediated_decay","non_stop_decay","protein_coding","ambiguous_orf","pseudogene","processed_pseudogene","polymorphic_pseudogene","transcribed_processed_pseudogene","transcribed_unprocessed_pseudogene","transcribed_unitary_pseudogene","translated_processed_pseudogene","translated_unprocessed_pseudogene","unitary_pseudogene","unprocessed_pseudogene","disrupted_domain"};
   int gt_size1 = 30;

   const char * gt_list2[100] = {"IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene","TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene","nonsense_mediated_decay","non_stop_decay","protein_coding","ambiguous_orf","disrupted_domain"};
   int gt_size2 = 14;

   int gt_size = gt_size2;
   const char ** gt_list;
   if (_used_pseudogene==1){
      gt_size = gt_size1;
      gt_list = gt_list1;
   }else{
      gt_size = gt_size2;
      gt_list = gt_list2;
   }

   std::string cur_g_name;
   std::string cur_g_type;
   std::string cur_g_str;
   while (std::getline(infile, line)){
      if (line.size()>0 && line[0]!='#'){
          m_substr_list = m_split_string(line, '\t');
          if (m_substr_list.size()>0){
             if (m_substr_list[0].size()<3 || m_substr_list[0].substr(0,3).compare("chr")!=0) {continue; }
             if (m_substr_list[0].find_first_of('_')!=std::string::npos){continue;}

             if (m_substr_list[2].compare("gene")==0){
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
                 cur_g_name.clear(); cur_g_type.clear();
                 while (std::getline(oss, cur_g_str, ' ')){
                     if (cur_g_str.compare("gene_name")==0){
                         std::getline(oss, cur_g_name, ' ');
                         cur_g_name = cur_g_name.substr(1, cur_g_name.size()-3);
                     }else if (cur_g_str.compare("gene_type")==0){
                         std::getline(oss, cur_g_type, ' ');
                         cur_g_type = cur_g_type.substr(1, cur_g_type.size()-3);
                     }
                 }
                 bool is_g=false;
                 for (int isg_i=0; isg_i<gt_size; isg_i++){
                    if (cur_g_type.compare(gt_list[isg_i])==0){
                       is_g = true;
                       break;
                    }
                 }
                 if (!is_g){continue;}
                 if (cur_g_name.size()==0){
                     std::cout<< "No gene_name find: " << line << std::endl;
                 }else{
                     gl_it = m_gene_list.find(cur_gr.chrn);
                     if (gl_it==m_gene_list.end()){
                         std::map<std::string, GenomicRegion> new_chr_map;
                         new_chr_map[cur_g_name] = cur_gr;
                         m_gene_list[cur_gr.chrn] = new_chr_map;
                     }else{
                         gl_it->second[cur_g_name] = cur_gr;
                     }
                 }
             }
          }
      }
   }
   infile.close();

   return m_gene_list;
}



