#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include <htslib/sam.h>

#include <memory>
#include <ctime>

#include <list>
#include <iostream>

#include <unordered_map>
#include <map>
#include <vector>

#include <sstream>
#include <fstream>

#include <climits>

#include <time.h>

#include <algorithm>    // std::find

#include "_gtf_struct_.h"
#include "_com_fun_.h"
#include "_com_structs_.h"


int m_check_gene_fusion(const char* in_bam_file, const char* in_gtf_file, const int min_len_ovlp, const int64_t _bin_size, const int _used_pseudogene, const int64_t _min_map_len, const int _used_secondary_alignment, const int _min_sup_read){
   std::map<std::string, std::map<std::string, std::shared_ptr<_gtf_entry_> > > m_gene_list; 
   std::map<std::string, std::string> _gid_to_gn;
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
   int retv = read_gtf(m_gene_list, _gid_to_gn, in_gtf_file, gt_list, gt_size);

   samFile * in_bam = NULL;
   bam_hdr_t * hdr;
   bam1_t * bam_one_alignment;
   bam1_core_t * bam_alignment_info;
   int sam_ret;

   uint32_t* m_cigar;
   int m_i_cigar; 
   int m_op;
   int m_len;

   int32_t ref_alg_pos;
   int32_t qry_alg_pos;
   int64_t num_del = 0;
   int64_t num_ins = 0;
   int64_t num_del_this = 0;
   int64_t num_ins_this = 0;

   uint64_t mapped_bases = 0;
   std::unordered_map<std::string, uint64_t> read_len_map;
   std::unordered_map<std::string, uint64_t> read_len_all;
   std::unordered_map<std::string, uint64_t>::iterator rl_map_it;
   std::list<uint32_t> read_len_list;
   std::list<uint32_t> read_len_map_list;
   std::list<uint32_t>::reverse_iterator rit;

   // key: query name; 
   std::unordered_map<std::string, std::vector<GenomicMapRegion> > multi_map_reads;
   std::unordered_map<std::string, std::vector<GenomicMapRegion> >::iterator multi_map_reads_it;

   fprintf(stdout, "Open sam file (%s).\n", in_bam_file);
   bam_one_alignment = bam_init1();
   bam_alignment_info = &bam_one_alignment->core;

   in_bam = sam_open(in_bam_file, "rb");
   if (NULL == in_bam){
      fprintf(stderr, "Error! Cannot open sam file (%s).\n", in_bam_file);
      abort();
   }

   hdr = sam_hdr_read(in_bam);
   std::string _this_q_name("");
   std::map<uint64_t, bool> _this_map_info;
   
   std::string m_qname;
   fprintf(stdout, "Read sam file (%s).\n", in_bam_file);
   fflush(stdout);
   clock_t _start_time = clock();
   while (sam_ret = sam_read1(in_bam, hdr, bam_one_alignment)>=0){
        m_qname = bam_get_qname(bam_one_alignment);
        //fprintf(stdout, "For read (%s).\n", m_qname.c_str());
        rl_map_it = read_len_all.find(m_qname);
        if (rl_map_it==read_len_all.end() || rl_map_it->second < bam_one_alignment->core.l_qseq){
            read_len_all[m_qname] = bam_one_alignment->core.l_qseq;
        } 

        if (bam_alignment_info->flag & BAM_FUNMAP) { continue; }
        if (_used_secondary_alignment==0 && (bam_alignment_info->flag & BAM_FSECONDARY)) { continue; }
        //if (_used_secondary_alignment==0 && (bam_alignment_info->flag & BAM_FSUPPLEMENTARY)) { continue; }

        rl_map_it = read_len_map.find(m_qname);
        if (rl_map_it==read_len_map.end() || rl_map_it->second < bam_one_alignment->core.l_qseq){
            read_len_map[m_qname] = bam_one_alignment->core.l_qseq;
        }

        //fprintf(stdout, "For read check (%s).\n", m_qname.c_str());
        if (_this_q_name.size()>0 && _this_q_name.compare(m_qname)!=0){
            num_del += num_del_this;
            num_ins += num_ins_this; 
            num_del_this = 0;
            num_ins_this = 0;
  
            if (multi_map_reads[_this_q_name].size()<2){
               multi_map_reads.erase(_this_q_name); 
            }
 
            if (read_len_map[_this_q_name] < _this_map_info.size()){
                fprintf(stderr, "Warning!! Mapped bases(%d) is more than total bases(%d) for %s.\n", read_len_map[_this_q_name], _this_map_info.size(), _this_q_name.c_str());
            }
            mapped_bases += _this_map_info.size();
            _this_map_info.clear();
        }
        _this_q_name = std::string(m_qname);
        
        multi_map_reads_it = multi_map_reads.find(m_qname);
        if (multi_map_reads_it==multi_map_reads.end()){
           multi_map_reads[m_qname] = std::vector<GenomicMapRegion>();
        }

        //fprintf(stdout, "For read pos (%s).\n", m_qname.c_str());
        ref_alg_pos = bam_alignment_info->pos;
        qry_alg_pos = 0;
        GenomicMapRegion _this_gmr;
        _this_gmr.chrn = hdr->target_name[bam_alignment_info->tid];
        _this_gmr.ref_start_pos = bam_alignment_info->pos;
        _this_gmr.qry_start_pos = 0;
        bool is_first_match = true;
        _this_gmr.ref_end_pos = 0;
        _this_gmr.qry_end_pos = 0;
        _this_gmr.qry_name = m_qname;
        if (bam_is_rev(bam_one_alignment) ) { _this_gmr.map_strand =1 ;}
        else { _this_gmr.map_strand = 0; }

        int64_t _this_read_len = 0;
        m_cigar = bam_get_cigar(bam_one_alignment);
        //fprintf(stdout, "For read pos_end (%s).\n", m_qname.c_str());
        for (m_i_cigar=0; m_i_cigar<bam_alignment_info->n_cigar; ++m_i_cigar){
           //fprintf(stdout, "For read cigar (%s).\n", m_qname.c_str());
           m_op = bam_cigar_op(m_cigar[m_i_cigar]);
           //fprintf(stdout, "For read cigar_len (%s).\n", m_qname.c_str());
           m_len = bam_cigar_oplen(m_cigar[m_i_cigar]);
           //fprintf(stdout, "\t%c%d\n", m_op, m_len);
           if (!(m_op==BAM_CDEL || m_op==BAM_CREF_SKIP || m_op==BAM_CPAD || m_op==BAM_CBACK)){
              _this_read_len += m_len;
           }
           switch (m_op){
              case BAM_CEQUAL:
              case BAM_CMATCH: // M
                   //fprintf(stdout, "%s%d", "M", m_len);
                   for (int ai=0; ai<m_len; ai++){
                      _this_map_info[qry_alg_pos+ai] = true;
                   }
                   if (is_first_match){
                      _this_gmr.qry_start_pos = qry_alg_pos;
                      is_first_match = false;
                   }
                   ref_alg_pos += m_len;
                   qry_alg_pos += m_len;
                   _this_gmr.ref_end_pos = ref_alg_pos;
                   _this_gmr.qry_end_pos = qry_alg_pos;
                   break;
              case BAM_CINS:  // I
                   //fprintf(stdout, "%s%d", "I", m_len);
                   num_ins_this += m_len;
                   /*for (int ai=0; ai<m_len; ai++){
                      _this_map_info[qry_alg_pos+ai] = true;
                   }*/
                   qry_alg_pos += m_len;
                   break;
              case BAM_CDEL:   // D
                   //fprintf(stdout, "%s%d", "D", m_len);
                   /*for (int ai=0; ai<m_len; ai++){
                      _this_map_info[ref_alg_pos+ai] = true;
                   }*/
                   num_del_this += m_len;
                   ref_alg_pos += m_len;
                   break;
              case BAM_CREF_SKIP:   // 
                   //fprintf(stdout, "%s%d", "R", m_len);
                   ref_alg_pos += m_len;
                   break;
              case BAM_CSOFT_CLIP:  //
                   //fprintf(stdout, "%s%d", "S", m_len);
                   //num_ins_this += m_len;
                   qry_alg_pos += m_len;
                   break;
              case BAM_CHARD_CLIP:  // 
                   //fprintf(stdout, "%s%d", "H", m_len);
                   //num_ins_this += m_len;
                   qry_alg_pos += m_len;
                   break;
              case BAM_CPAD:  //
                   //fprintf(stdout, "%s%d", "P", m_len);
                   break;
              case BAM_CDIFF: //
                   //fprintf(stdout, "%s%d", "X", m_len);
                   for (int ai=0; ai<m_len; ai++){
                      _this_map_info[qry_alg_pos+ai] = true;
                   }
                   if (is_first_match){
                      _this_gmr.qry_start_pos = qry_alg_pos;
                      is_first_match = false;
                   }
                   ref_alg_pos += m_len;
                   qry_alg_pos += m_len;
                   _this_gmr.ref_end_pos = ref_alg_pos;
                   _this_gmr.qry_end_pos = qry_alg_pos;
                   break;
              default:
                   fprintf(stderr, "Unknow cigar %d:%d\n", m_op, m_len);              
           }
        }
        if (_this_gmr.map_strand == 1){
           int64_t _rd_ind = _this_read_len - _this_gmr.qry_start_pos -1;
           _this_gmr.qry_start_pos = _this_read_len - _this_gmr.qry_end_pos -1;
           _this_gmr.qry_end_pos = _rd_ind;
        }
        //if (_this_gmr.qry_name.compare("35e82b56-9386-4144-94a3-15e560f45f40")==0 || _this_gmr.qry_name.compare("90bf2a20-9476-4c26-9a86-bb69bfa6a980")==0 || _this_gmr.qry_name.compare("a61c8185-46f7-4d00-8c1a-69a1b3b6efdb")==0 ){
        //   std::cout<<" Test "<<_this_gmr.qry_name<<":"<<_this_gmr.qry_start_pos<<"-"<<_this_gmr.qry_end_pos<<" " <<_this_gmr.map_strand<<"/"<<_this_gmr.chrn<<":"<<_this_gmr.ref_start_pos<<"-"<<_this_gmr.ref_end_pos<<std::endl;
        //}

        if (_min_map_len<_this_gmr.qry_end_pos-_this_gmr.qry_start_pos){
           multi_map_reads[m_qname].push_back(_this_gmr);
        }
   }
   fprintf(stdout, "Time consumed to read bam: %.2fs\n", (double)(clock() - _start_time)/CLOCKS_PER_SEC);

   if (_this_q_name.size()>0 && _this_q_name.compare(m_qname)!=0){
      num_del += num_del_this;
      num_ins += num_ins_this;
      num_del_this = 0;
      num_ins_this = 0;

      if (multi_map_reads[_this_q_name].size()<2){
          multi_map_reads.erase(_this_q_name);
      }

      if (read_len_map[_this_q_name] < _this_map_info.size()){
          fprintf(stderr, "Warning!! Mapped bases(%d) is more than total bases(%d) for %s.\n", read_len_map[_this_q_name], _this_map_info.size(), _this_q_name.c_str());
      }
      mapped_bases += _this_map_info.size();
      _this_map_info.clear();
   }
   _this_q_name = std::string(m_qname);

   uint64_t tot_all_bases = 0;
   uint64_t tot_map_bases = 0;
   read_len_map_list.clear();
   read_len_list.clear();
   for (rl_map_it=read_len_map.begin(); rl_map_it!=read_len_map.end(); rl_map_it++){
       read_len_map_list.push_back(rl_map_it->second);
       tot_map_bases += rl_map_it->second;
   }
   for (rl_map_it=read_len_all.begin(); rl_map_it!=read_len_all.end(); rl_map_it++){
       read_len_list.push_back(rl_map_it->second);
       tot_all_bases += rl_map_it->second;
   }

   uint64_t total_reads = read_len_all.size();
   uint64_t maped_reads = read_len_map.size();   

   read_len_map_list.sort();
   read_len_list.sort();  

   std::map<uint64_t, uint64_t> _all_np;
   std::map<uint64_t, uint64_t> _map_np;
   for(int _np=10; _np<100; _np+=10){
      _all_np[_np] = 0;
      _map_np[_np] = 0;
   }

   uint64_t csize = 0; 
   for (rit=read_len_map_list.rbegin(); rit!=read_len_map_list.rend(); ++rit){
      csize += *rit;
      for(int _np=10; _np<100; _np+=10){
         if (_map_np[_np]==0 && ((double(csize)/tot_map_bases)>=_np/float(100))){
            _map_np[_np] = *rit;
         }
      }
   }

   csize = 0;
   for (rit=read_len_list.rbegin(); rit!=read_len_list.rend(); ++rit){
      csize += *rit;
      for(int _np=10; _np<100; _np+=10){
         if (_all_np[_np]==0 && ((double(csize)/tot_all_bases)>=_np/float(100))){
            _all_np[_np] = *rit;
         }
      }
   }

   fprintf(stdout, "For all reads: reads=%llu/%llu(%.2f) bases=%llu/%llu(%.2f) ins=%llu(%.2f) del=%llu(%.2f)\n", 
                   read_len_map_list.size(),read_len_list.size(),((double)(read_len_map_list.size())/read_len_list.size())*100,
                   mapped_bases,tot_all_bases,((double)mapped_bases/tot_all_bases)*100,
                   num_ins,(double)num_ins*100/tot_all_bases, num_del,(double)num_del*100/tot_all_bases );
   for(int _np=10; _np<100; _np+=10){
       fprintf(stdout, "\t N%d: %d\n", _np, _all_np[_np]);
   }
   fprintf(stdout, "For mapped reads: bases=%llu/%llu(%.2f) ins=%llu(%.2f) del=%llu(%.2f)\n", 
                   mapped_bases,tot_map_bases,((double)mapped_bases/tot_map_bases)*100,
                   num_ins,(double)num_ins*100/tot_map_bases, num_del,(double)num_del*100/tot_map_bases );
   for(int _np=10; _np<100; _np+=10){
       fprintf(stdout, "\t N%d: %d\n", _np, _map_np[_np]);
   }
   fprintf(stdout, "\n");

   // for gene fusion
   fprintf(stdout, "# reads with multipe mapping = %llu/%llu(%.2f)\n", multi_map_reads.size(), read_len_list.size(), ((double)(multi_map_reads.size())/read_len_list.size())*100);
   fflush(stdout);
   // key: gene pair
   std::map<std::string, GeneFusionCand> _fg_cand;
   std::map<std::string, GeneFusionCand>::iterator _fg_it;
   std::map<std::string, std::vector<GenomicMapRegion> >::iterator _g_map_it;
   std::map<std::string, int64_t> _gene_in_multi_map;
   clock_t _start_t = std::clock();
   int64_t _h_i = 0;
   for (multi_map_reads_it=multi_map_reads.begin(); multi_map_reads_it!=multi_map_reads.end(); multi_map_reads_it++){
      _h_i += 1;
      if (_h_i%50000==0 || _h_i<5){
          std::cout<<"Handle ="<<_h_i<<"/"<<multi_map_reads.size()<<" "<<_fg_cand.size()<<" Time-elapsed="<<int(double(std::clock()-_start_t)*100/CLOCKS_PER_SEC)/100.00<<std::endl;
      }
      //key: gene name
      std::map<std::string, std::vector<GenomicMapRegion> > _g_list;
      // get all alignment with a gene
      for(std::vector<GenomicMapRegion>::iterator _gmr_it=multi_map_reads_it->second.begin(); _gmr_it!=multi_map_reads_it->second.end(); _gmr_it++){
         std::map<std::string, std::map<std::string, std::shared_ptr<_gtf_entry_> >::iterator > _ovlp_rec;
         std::map<std::string, int64_t > _ovlp_len;
         for (std::map<std::string, std::map<std::string, std::shared_ptr<_gtf_entry_> > >::iterator gl_it=m_gene_list.begin(); gl_it!=m_gene_list.end(); gl_it++){
            if(gl_it->first.compare(_gmr_it->chrn)!=0) { continue; }
            for(std::map<std::string, std::shared_ptr<_gtf_entry_> >::iterator chr_g_it=gl_it->second.begin(); chr_g_it!=gl_it->second.end(); chr_g_it++){
                //if ((_gmr_it->qry_name.compare("35e82b56-9386-4144-94a3-15e560f45f40")==0 || _gmr_it->qry_name.compare("90bf2a20-9476-4c26-9a86-bb69bfa6a980")==0 || _gmr_it->qry_name.compare("a61c8185-46f7-4d00-8c1a-69a1b3b6efdb")==0) && (chr_g_it->second->get_name().compare("PSMA6")==0 || chr_g_it->second->get_name().compare("VMP1")==0)){
                //   std::cout<<" Test "<<_gmr_it->qry_name<<":"<<_gmr_it->qry_start_pos<<"-"<<_gmr_it->qry_end_pos<<" " <<_gmr_it->map_strand<<"/"<<_gmr_it->chrn<<":"<<_gmr_it->ref_start_pos<<"-"<<_gmr_it->ref_end_pos<<" "<<chr_g_it->second->get_ovlp(*_gmr_it)<<"/"<<chr_g_it->second->get_coding_ovlp(*_gmr_it)<<" "<<chr_g_it->second->_to_string()<<std::endl;
                //}

                if (chr_g_it->second->get_ovlp(*_gmr_it) > min_len_ovlp && chr_g_it->second->get_coding_ovlp(*_gmr_it) > min_len_ovlp){
                    if (_ovlp_rec.find(chr_g_it->first)==_ovlp_rec.end()){
                       _ovlp_rec[chr_g_it->first] = chr_g_it;
                       _ovlp_len[chr_g_it->first] = chr_g_it->second->get_coding_ovlp(*_gmr_it);
                    }else{
                       std::cout<<"Warning!!! Duplicate genes " << chr_g_it->first <<std::endl;
                    }
                    continue;

                    _g_map_it = _g_list.find(chr_g_it->first);
                    if (_g_map_it==_g_list.end()){
                       _g_list[chr_g_it->first] = std::vector<GenomicMapRegion>();
                    }
                    _g_list[chr_g_it->first].push_back(*_gmr_it);
                    if(_gene_in_multi_map.find(chr_g_it->first)==_gene_in_multi_map.end()){
                       _gene_in_multi_map[chr_g_it->first] = 0;
                    }
                    _gene_in_multi_map[chr_g_it->first] += 1;
                }
            }
         }
         int64_t _max_l = 0;
         std::string _max_l_k;
         for(std::map<std::string, int64_t >::iterator _ovlp_l_i=_ovlp_len.begin(); _ovlp_l_i!=_ovlp_len.end(); _ovlp_l_i++){
            if (_ovlp_l_i->second>_max_l){
               _max_l = _ovlp_l_i->second;
               _max_l_k = _ovlp_l_i->first;
            }
         }
         if (_max_l>min_len_ovlp){
            std::map<std::string, std::shared_ptr<_gtf_entry_> >::iterator chr_g_it= _ovlp_rec[_max_l_k];
            _g_map_it = _g_list.find(chr_g_it->first);
            if (_g_map_it==_g_list.end()){
                _g_list[chr_g_it->first] = std::vector<GenomicMapRegion>();
            }
            _g_list[chr_g_it->first].push_back(*_gmr_it);
            if(_gene_in_multi_map.find(chr_g_it->first)==_gene_in_multi_map.end()){
                _gene_in_multi_map[chr_g_it->first] = 0;
            }
            _gene_in_multi_map[chr_g_it->first] += 1;
         }
      }
      // less than 2 genes; ignore.
      if (_g_list.size()<2){ continue; }
      //if (_g_list.size()<_min_sup_read){ continue; }

      std::map<std::string, std::vector<GenomicMapRegion> >::iterator g_l_it1 = _g_list.begin();   
      std::map<std::string, std::vector<GenomicMapRegion> >::iterator g_l_it2 = _g_list.begin();
      std::map<std::string, std::vector<GenomicMapRegion> >::iterator g_l_it_t1; 
      std::map<std::string, std::vector<GenomicMapRegion> >::iterator g_l_it_t2;
      g_l_it2++;
      for(; g_l_it2!=_g_list.end(); g_l_it1++,g_l_it2++){
          // check gene overlap; 
          // ignore those gene pair whose two genes are overlapped.
          /*if (m_gene_list[g_l_it1->second[0].chrn][g_l_it1->first].chrn.compare(m_gene_list[g_l_it2->second[0].chrn][g_l_it2->first].chrn)==0){
             int64_t max_start =  m_gene_list[g_l_it1->second[0].chrn][g_l_it1->first].start_pos > m_gene_list[g_l_it2->second[0].chrn][g_l_it2->first].start_pos ? m_gene_list[g_l_it1->second[0].chrn][g_l_it1->first].start_pos : m_gene_list[g_l_it2->second[0].chrn][g_l_it2->first].start_pos;
             int64_t min_end = m_gene_list[g_l_it1->second[0].chrn][g_l_it1->first].end_pos < m_gene_list[g_l_it2->second[0].chrn][g_l_it2->first].end_pos ? m_gene_list[g_l_it1->second[0].chrn][g_l_it1->first].end_pos : m_gene_list[g_l_it2->second[0].chrn][g_l_it2->first].end_pos;
             if (min_end - max_start > 10){
                 continue;
             }else{
                if (m_gene_list[g_l_it1->second[0].chrn][g_l_it1->first].get_ovlp(m_gene_list[g_l_it2->second[0].chrn][g_l_it2->first]) > 10){
                   std::cout<<"\t gene ovlp"<<m_gene_list[g_l_it1->second[0].chrn][g_l_it1->first]._to_string()<<" "<<m_gene_list[g_l_it2->second[0].chrn][g_l_it2->first]._to_string()<< " ovlp="<<(m_gene_list[g_l_it1->second[0].chrn][g_l_it1->first].get_ovlp(m_gene_list[g_l_it2->second[0].chrn][g_l_it2->first])<<std::endl;
                }
             }
          }*/
          //std::cout<<" db "<<g_l_it1->second[0].chrn<<":"<<g_l_it1->first<< " " << g_l_it2->second[0].chrn<<":"<<g_l_it2->first << " "<<m_gene_list[g_l_it1->second[0].chrn][g_l_it1->first]->get_ovlp(m_gene_list[g_l_it2->second[0].chrn][g_l_it2->first]) <<std::endl;
          if (m_gene_list[g_l_it1->second[0].chrn][g_l_it1->first]->get_ovlp(m_gene_list[g_l_it2->second[0].chrn][g_l_it2->first]) > 10){
             continue;
          }

          std::string _gene_pair;
          if (g_l_it1->first>g_l_it2->first){
             _gene_pair = g_l_it2->first + g_l_it1->first;
             g_l_it_t1 = g_l_it2; 
             g_l_it_t2 = g_l_it1;
          }else if (g_l_it1->first<g_l_it2->first){
             _gene_pair = g_l_it1->first + g_l_it2->first;
             g_l_it_t1 = g_l_it1;
             g_l_it_t2 = g_l_it2;
          }else{
             fprintf(stderr, "Warning!!! Two genes are same: %s vs %s\n", g_l_it1->first.c_str(), g_l_it2->first.c_str());
             continue;
          }
          _fg_it = _fg_cand.find(_gene_pair);
          if (_fg_it==_fg_cand.end()){
              GeneFusionCand _new_fg;
              if (g_l_it1->first>g_l_it2->first){
                  _new_fg.g1 = g_l_it2->first;
                  _new_fg.g2 = g_l_it1->first;
              }else{
                  _new_fg.g1 = g_l_it1->first;
                  _new_fg.g2 = g_l_it2->first;
              }
              _fg_cand[_gene_pair] = _new_fg;
          }
         
          std::map<std::string, std::vector<GenomicMapRegion> >::iterator _q_m_it;
          // cand*: key------query name
          _q_m_it = _fg_cand[_gene_pair].cand1.find(multi_map_reads_it->first); 
          if (_q_m_it==_fg_cand[_gene_pair].cand1.end()){
             _fg_cand[_gene_pair].cand1[multi_map_reads_it->first] = std::vector<GenomicMapRegion>();
          }
          for(std::vector<GenomicMapRegion>::iterator _gp_gmr_it=g_l_it_t1->second.begin(); _gp_gmr_it!=g_l_it_t1->second.end(); _gp_gmr_it++){
             _fg_cand[_gene_pair].cand1[multi_map_reads_it->first].push_back(*_gp_gmr_it);
          }

          _q_m_it = _fg_cand[_gene_pair].cand2.find(multi_map_reads_it->first);
          if (_q_m_it==_fg_cand[_gene_pair].cand2.end()){
             _fg_cand[_gene_pair].cand2[multi_map_reads_it->first] = std::vector<GenomicMapRegion>();
          }
          for(std::vector<GenomicMapRegion>::iterator _gp_gmr_it=g_l_it_t2->second.begin(); _gp_gmr_it!=g_l_it_t2->second.end(); _gp_gmr_it++){
             _fg_cand[_gene_pair].cand2[multi_map_reads_it->first].push_back(*_gp_gmr_it);
          }
      }
   }

   std::cout<<"Finish _fg_cand[_gene_pair="<<_fg_cand.size()<<std::endl;
   fflush(stdout);

   std::vector<NumStrPair> _gp_rank;
   std::map<std::string, std::string> _gp_pair_map;
   for (_fg_it = _fg_cand.begin(); _fg_it != _fg_cand.end(); _fg_it++){
      int64_t _g_str_cand1 = LLONG_MAX;
      int64_t _g_str_cand2 = LLONG_MAX;
      int64_t _g_end_cand1 = 0;
      int64_t _g_end_cand2 = 0;
      
      for(int _c_i = 1; _c_i<3; _c_i++){
         std::map<std::string, std::vector<GenomicMapRegion> > & _thi_cand = (_c_i==1?_fg_it->second.cand1:_fg_it->second.cand2);
         std::map<std::string, std::vector<GenomicMapRegion> >::iterator _qry_it = _thi_cand.begin();
         for (; _qry_it!=_thi_cand.end(); _qry_it++){
            std::vector<GenomicMapRegion>::iterator _q_grm_it = _qry_it->second.begin();
            for (; _q_grm_it!= _qry_it->second.end(); _q_grm_it++){
               if (_c_i==1){
                  if (_g_str_cand1>_q_grm_it->ref_start_pos){ _g_str_cand1 = _q_grm_it->ref_start_pos; }
                  if (_g_end_cand1<_q_grm_it->ref_end_pos){ _g_end_cand1 = _q_grm_it->ref_end_pos; }
               }else{
                  if (_g_str_cand2>_q_grm_it->ref_start_pos){ _g_str_cand2 = _q_grm_it->ref_start_pos; }
                  if (_g_end_cand2<_q_grm_it->ref_end_pos){ _g_end_cand2 = _q_grm_it->ref_end_pos; }
               }
            }
         }
      }
     
      // key: pos_bin 
      std::map<std::string, GeneFusionCand > _d_cand_pair;
      if (true){
         std::map<int, std::map<std::string, int> >::iterator _g_p_sup_reads_it;
         std::map<std::string, int>::iterator _g_p_sup_it;

         std::map<std::string, std::vector<GenomicMapRegion> > & _thi_cand = _fg_it->second.cand1;
         std::map<std::string, std::vector<GenomicMapRegion> > & _oth_cand = _fg_it->second.cand2;
         std::map<std::string, std::vector<GenomicMapRegion> >::iterator _qry_it = _thi_cand.begin();
         for (; _qry_it!=_thi_cand.end(); _qry_it++){
            std::vector<GenomicMapRegion>::iterator _q_grm_it = _qry_it->second.begin();
            for (; _q_grm_it!= _qry_it->second.end(); _q_grm_it++){
               std::vector<GenomicMapRegion>::iterator _oq_grm_it = _oth_cand[_qry_it->first].begin();
               for(; _oq_grm_it!= _oth_cand[_qry_it->first].end(); _oq_grm_it++){
                  int64_t max_str = _q_grm_it->qry_start_pos > _oq_grm_it->qry_start_pos ? _q_grm_it->qry_start_pos : _oq_grm_it->qry_start_pos;
                  int64_t min_end = _q_grm_it->qry_end_pos   > _oq_grm_it->qry_end_pos   ? _oq_grm_it->qry_end_pos  : _q_grm_it->qry_end_pos;
                  if (min_end - max_str > 50 || min_end - max_str < -20){ continue; }
                  if ((min_end - max_str)/double(_q_grm_it->qry_end_pos - _q_grm_it->qry_start_pos)>0.4){ continue; }
                  if ((min_end - max_str)/double(_oq_grm_it->qry_end_pos - _oq_grm_it->qry_start_pos)>0.4){ continue; }

                  /*if (_oq_grm_it->qry_name.compare("0f1b38a9-96b2-4afc-9544-958712bebd00")==0 || _oq_grm_it->qry_name.compare("642a243d-d6ac-4db0-8a44-3aa3a01dc26c")==0 || _oq_grm_it->qry_name.compare("716353fb-78ac-4cf2-80be-5a5044869374")==0) {
                      std::cout<<"Test.1: "<<(_q_grm_it->map_strand==0?"+":"-")<<_q_grm_it->chrn<<":"<<_q_grm_it->ref_start_pos<<"-"<<_q_grm_it->ref_end_pos<<"/"<<_q_grm_it->qry_name<<":"<<_q_grm_it->qry_start_pos<<"-"<<_q_grm_it->qry_end_pos<<" >>> "<<max_str<<" "<<min_end<<std::endl;
                      std::cout<<"Test.2: "<<(_oq_grm_it->map_strand==0?"+":"-")<<_oq_grm_it->chrn<<":"<<_oq_grm_it->ref_start_pos<<"-"<<_oq_grm_it->ref_end_pos<<"/"<<_oq_grm_it->qry_name<<":"<<_oq_grm_it->qry_start_pos<<"-"<<_oq_grm_it->qry_end_pos<<" >>> "<<min_end-max_str<<std::endl;
                  }*/

                  if (_q_grm_it->qry_end_pos < _q_grm_it->qry_start_pos) {
                      std::cout<<"Error"<<(_q_grm_it->map_strand==0?"+":"-")<<_q_grm_it->chrn<<":"<<_q_grm_it->ref_start_pos<<"-"<<_q_grm_it->ref_end_pos<<"/"<<_q_grm_it->qry_name<<":"<<_q_grm_it->qry_start_pos<<"-"<<_q_grm_it->qry_end_pos<<std::endl;
                  }
                  if (_oq_grm_it->qry_end_pos < _oq_grm_it->qry_start_pos){
                      std::cout<<"Error"<<(_oq_grm_it->map_strand==0?"+":"-")<<_oq_grm_it->chrn<<":"<<_oq_grm_it->ref_start_pos<<"-"<<_oq_grm_it->ref_end_pos<<"/"<<_oq_grm_it->qry_name<<":"<<_oq_grm_it->qry_start_pos<<"-"<<_oq_grm_it->qry_end_pos<<std::endl;
                  }

                  int64_t pos_i_1 = -1;
                  int64_t pos_i_2 = -1;
                  int64_t bk1 = -1;
                  int64_t bk2 = -1;
                  if (_q_grm_it->map_strand==0 && _oq_grm_it->map_strand==0){
                     if (_q_grm_it->qry_start_pos < _oq_grm_it->qry_start_pos){
                        pos_i_1 = (_q_grm_it->ref_end_pos - _g_str_cand1)/(_bin_size/2); 
                        pos_i_2 = ( _oq_grm_it->ref_start_pos - _g_str_cand2)/(_bin_size/2);
                        bk1 = _q_grm_it->ref_end_pos;
                        bk2 = _oq_grm_it->ref_start_pos;
                     }else{
                        pos_i_1 = (_q_grm_it->ref_start_pos - _g_str_cand1)/(_bin_size/2);
                        pos_i_2 = ( _oq_grm_it->ref_end_pos - _g_str_cand2)/(_bin_size/2);
                        bk1 = _q_grm_it->ref_start_pos;
                        bk2 = _oq_grm_it->ref_end_pos;
                     }
                  }else if (_q_grm_it->map_strand==0 && _oq_grm_it->map_strand==1){
                     if (_q_grm_it->qry_start_pos < _oq_grm_it->qry_start_pos){
                        pos_i_1 = (_q_grm_it->ref_end_pos - _g_str_cand1)/(_bin_size/2);
                        pos_i_2 = ( _oq_grm_it->ref_end_pos - _g_str_cand2)/(_bin_size/2);
                        bk1 = _q_grm_it->ref_end_pos;
                        bk2 =  _oq_grm_it->ref_end_pos;
                     }else{
                        pos_i_1 = (_q_grm_it->ref_start_pos - _g_str_cand1)/(_bin_size/2);
                        pos_i_2 = ( _oq_grm_it->ref_start_pos - _g_str_cand2)/(_bin_size/2);
                        bk1 = _q_grm_it->ref_start_pos;
                        bk2 = _oq_grm_it->ref_start_pos;
                     }
                  }else if (_q_grm_it->map_strand==1 && _oq_grm_it->map_strand==0){
                     if (_q_grm_it->qry_start_pos < _oq_grm_it->qry_start_pos){
                        pos_i_1 = (_q_grm_it->ref_start_pos - _g_str_cand1)/(_bin_size/2);
                        pos_i_2 = ( _oq_grm_it->ref_start_pos - _g_str_cand2)/(_bin_size/2);
                        bk1 = _q_grm_it->ref_start_pos;
                        bk2 = _oq_grm_it->ref_start_pos;
                     }else{
                        pos_i_1 = (_q_grm_it->ref_end_pos - _g_str_cand1)/(_bin_size/2);
                        pos_i_2 = ( _oq_grm_it->ref_end_pos - _g_str_cand2)/(_bin_size/2);
                        bk1 = _q_grm_it->ref_end_pos;
                        bk2 =  _oq_grm_it->ref_end_pos;
                     }
                  }else if (_q_grm_it->map_strand==1 && _oq_grm_it->map_strand==1){
                     if (_q_grm_it->qry_start_pos < _oq_grm_it->qry_start_pos){
                        pos_i_1 = (_q_grm_it->ref_start_pos - _g_str_cand1)/(_bin_size/2);
                        pos_i_2 = ( _oq_grm_it->ref_end_pos - _g_str_cand2)/(_bin_size/2);
                        bk1 = _q_grm_it->ref_start_pos;
                        bk2 =  _oq_grm_it->ref_end_pos;
                     }else{
                        pos_i_1 = (_q_grm_it->ref_end_pos - _g_str_cand1)/(_bin_size/2);
                        pos_i_2 = ( _oq_grm_it->ref_start_pos - _g_str_cand2)/(_bin_size/2);
                        bk1 = _q_grm_it->ref_end_pos;
                        bk2 =  _oq_grm_it->ref_start_pos;
                     }
                  }else{
                      fprintf(stderr, "Warning!!! Wrong strand: %d vs %d for %s\n", _q_grm_it->map_strand, _oq_grm_it->map_strand, _qry_it->first.c_str());
                      continue;
                  }

                  if (bk1 < m_gene_list[_q_grm_it->chrn][_fg_it->second.g1]->get_gr_start_pos() || bk1 > m_gene_list[_q_grm_it->chrn][_fg_it->second.g1]->get_gr_end_pos() ) {
                     if (m_gene_list[_q_grm_it->chrn][_fg_it->second.g1]->get_name().compare("PSMA6")==0 || m_gene_list[_q_grm_it->chrn][_fg_it->second.g1]->get_name().compare("VMP1")==0){
                        std::cout<<"Test2 "<<m_gene_list[_q_grm_it->chrn][_fg_it->second.g1]->_to_string()<<" >>> "<<bk1<<std::endl;
                     }
                     continue;
                  }
                  if (bk2 < m_gene_list[_oq_grm_it->chrn][_fg_it->second.g2]->get_gr_start_pos() || bk2 > m_gene_list[_oq_grm_it->chrn][_fg_it->second.g2]->get_gr_end_pos() ) {
                     if (m_gene_list[_oq_grm_it->chrn][_fg_it->second.g2]->get_name().compare("PSMA6")==0 || m_gene_list[_oq_grm_it->chrn][_fg_it->second.g2]->get_name().compare("VMP1")==0){
                        std::cout<<"Test2 "<<m_gene_list[_oq_grm_it->chrn][_fg_it->second.g2]->_to_string()<<" >>> "<<bk2 <<std::endl;
                     }
                     continue;
                  }

                  _q_grm_it->other_pos = bk1;
                  _oq_grm_it->other_pos = bk2;

                  for (int64_t _add_i_1 = pos_i_1-1; _add_i_1<pos_i_1+1; _add_i_1++){
                      if (_add_i_1<-1){
                         fprintf(stderr, "Warning!!! Wrong1 %d:%llu-%llu %d:%llu-%llu %llu %llu", _q_grm_it->map_strand, _q_grm_it->ref_start_pos, _q_grm_it->ref_end_pos, _oq_grm_it->map_strand, _oq_grm_it->ref_start_pos, _oq_grm_it->ref_end_pos, _q_grm_it->qry_start_pos, _oq_grm_it->qry_start_pos); continue;
                      }else if (_add_i_1<0){ continue; }
                      for (int64_t _add_i_2 = pos_i_2-1; _add_i_2<pos_i_2+1; _add_i_2++){
                          if (_add_i_2<-1){
                             fprintf(stderr, "Warning!!! Wrong2 %d:%llu-%llu %d:%llu-%llu %llu %llu\n", _q_grm_it->map_strand, _q_grm_it->ref_start_pos, _q_grm_it->ref_end_pos, _oq_grm_it->map_strand, _oq_grm_it->ref_start_pos, _oq_grm_it->ref_end_pos, _q_grm_it->qry_start_pos, _oq_grm_it->qry_start_pos); continue;
                          }else if (_add_i_2<0){ continue; }
                      
                         std::string _pos_pair = std::string(std::to_string(_add_i_1))+"_"+std::string(std::to_string(_add_i_2));
                         if (_d_cand_pair.find(_pos_pair) == _d_cand_pair.end()){
                             GeneFusionCand _p_gfc;
                             _p_gfc.g1 = _fg_it->second.g1;
                             _p_gfc.g2 = _fg_it->second.g2;
                             _d_cand_pair[_pos_pair] = _p_gfc;
                         }
                         if (_d_cand_pair[_pos_pair].cand1.find(_qry_it->first) == _d_cand_pair[_pos_pair].cand1.end()){
                             _d_cand_pair[_pos_pair].cand1[_qry_it->first] = std::vector<GenomicMapRegion>();
                         }
                         _d_cand_pair[_pos_pair].cand1[_qry_it->first].push_back(*_q_grm_it);
                         if (_d_cand_pair[_pos_pair].cand2.find(_qry_it->first) == _d_cand_pair[_pos_pair].cand2.end()){
                             _d_cand_pair[_pos_pair].cand2[_qry_it->first] = std::vector<GenomicMapRegion>();
                         }
                         _d_cand_pair[_pos_pair].cand2[_qry_it->first].push_back(*_oq_grm_it);
                      }
                  }
               }
            } 
         }
      }
      // 
      int max_sup_1 = 0;
      int max_sup_2 = 0;
      int64_t tot_sup_1=0, tot_sup_2=0;
      std::map<std::string, GeneFusionCand >::iterator _d_dp_it;
      for (std::map<std::string, GeneFusionCand >::iterator _d_dp_it_go=_d_cand_pair.begin(); _d_dp_it_go!=_d_cand_pair.end(); _d_dp_it_go++){
          tot_sup_1 += _d_dp_it_go->second.cand1.size();  tot_sup_2 += _d_dp_it_go->second.cand2.size();
          if (_d_dp_it_go->second.cand1.size() > max_sup_1){ max_sup_1 = _d_dp_it_go->second.cand1.size(); _d_dp_it=_d_dp_it_go; }
          if (_d_dp_it_go->second.cand2.size() > max_sup_2){ max_sup_2 = _d_dp_it_go->second.cand2.size(); _d_dp_it=_d_dp_it_go; } 
      }
      if (max_sup_1>_min_sup_read-1 && max_sup_2>_min_sup_read-1){
         if (_d_dp_it->second.cand1.size()>_min_sup_read-1 || _d_dp_it->second.cand2.size()>_min_sup_read-1){
            std::ostringstream oss_tostr;
            oss_tostr<<"GF\t"<<_gid_to_gn[_fg_it->second.g1]<<":"<<_gid_to_gn[_fg_it->second.g2]<<" "<<max_sup_1<<" "<<max_sup_2<<" supporting reads="<<_d_dp_it->second.cand1.size()<<"/"<<_d_dp_it->second.cand2.size(); //_com_qry_name.size(); 
            std::vector<int64_t> pos_1, pos_2;
            int64_t pos_1_left=0, pos_1_right=0;
            int64_t pos_2_left=0, pos_2_right=0;
            for(std::map<std::string, std::vector<GenomicMapRegion> >::iterator _i_gmr_=_d_dp_it->second.cand1.begin(); _i_gmr_!=_d_dp_it->second.cand1.end(); _i_gmr_++){
               for (std::vector<GenomicMapRegion>::iterator _i=_i_gmr_->second.begin(); _i!=_i_gmr_->second.end(); _i++){
                  pos_1.push_back(_i->other_pos);
                  if (_i->other_pos == _i->ref_start_pos){ pos_1_left += 1; }
                  else if (_i->other_pos == _i->ref_end_pos) { pos_1_right += 1; }
                  else {std::cout<<"Error1 not start("<<_i->ref_start_pos<<") or end("<<_i->ref_end_pos<<") "<<_i->other_pos<<std::endl;}
               }
               for (std::vector<GenomicMapRegion>::iterator _i=_d_dp_it->second.cand2[_i_gmr_->first].begin(); _i!=_d_dp_it->second.cand2[_i_gmr_->first].end(); _i++){
                  pos_2.push_back(_i->other_pos);
                  if (_i->other_pos == _i->ref_start_pos){ pos_2_left += 1; }
                  else if (_i->other_pos == _i->ref_end_pos) { pos_2_right += 1; }
                  else {std::cout<<"Error2 not start("<<_i->ref_start_pos<<") or end("<<_i->ref_end_pos<<") "<<_i->other_pos<<std::endl;}
               }
            }
            std::map<std::string, std::vector<GenomicMapRegion> >::iterator _i_1=_d_dp_it->second.cand1.begin();
            std::map<std::string, std::vector<GenomicMapRegion> >::iterator _i_2=_d_dp_it->second.cand2.begin();
            oss_tostr<<" "<<_i_1->second[0].chrn<<":"<<int64_t(std::accumulate(std::begin(pos_1), std::end(pos_1), 0.0) / pos_1.size())<<" "<<pos_1_left<<"/"<<pos_1_right<<":"<<tot_sup_1<<":"<<_gene_in_multi_map[_fg_it->second.g1]<<" "<<_i_2->second[0].chrn<<":"<<int64_t(std::accumulate(std::begin(pos_2), std::end(pos_2), 0.0) / pos_2.size())<<" "<<pos_2_left<<"/"<<pos_2_right<<":"<<tot_sup_2<<":"<<_gene_in_multi_map[_fg_it->second.g2]<<"\n";
            for(std::map<std::string, std::vector<GenomicMapRegion> >::iterator _i=_d_dp_it->second.cand1.begin(); _i!=_d_dp_it->second.cand1.end(); _i++){
               oss_tostr<<"\t"<<_i->second[0].other_pos<<"("<<(_i->second[0].map_strand==0?"+":"-")<<_i->second[0].chrn<<":"<<_i->second[0].ref_start_pos<<"-"<<_i->second[0].ref_end_pos<<"/"<<_i->second[0].qry_name<<":"<<_i->second[0].qry_start_pos<<"-"<<_i->second[0].qry_end_pos<<")"<<_i->second.size()<<" ";
               oss_tostr<<_d_dp_it->second.cand2[_i->first][0].other_pos<<"("<<(_d_dp_it->second.cand2[_i->first][0].map_strand==0?"+":"-")<<_d_dp_it->second.cand2[_i->first][0].chrn<<":"<<_d_dp_it->second.cand2[_i->first][0].ref_start_pos<<"-"<<_d_dp_it->second.cand2[_i->first][0].ref_end_pos<<"/"<<_d_dp_it->second.cand2[_i->first][0].qry_start_pos<<"-"<<_d_dp_it->second.cand2[_i->first][0].qry_end_pos<<")"<<_d_dp_it->second.cand2[_i->first].size();
               oss_tostr<<"\n";  
            }
            
            for(; _i_1!=_d_dp_it->second.cand1.end(); _i_1++){
               for(std::vector<GenomicMapRegion>::iterator _i=_i_1->second.begin(); _i!=_i_1->second.end(); _i++){
                  std::cout<<"\t 1.INFO"<<_i->other_pos<<"("<<(_i->map_strand==0?"+":"-")<<_i->chrn<<":"<<_i->ref_start_pos<<"-"<<_i->ref_end_pos<<"/"<<_i->qry_name<<":"<<_i->qry_start_pos<<"-"<<_i->qry_end_pos<<") "<<std::endl;
               }
            }
            for(; _i_2!=_d_dp_it->second.cand2.end(); _i_2++){
               for(std::vector<GenomicMapRegion>::iterator _i=_i_2->second.begin(); _i!=_i_2->second.end(); _i++){
                  std::cout<<"\t 2.INFO"<<_i->other_pos<<"("<<(_i->map_strand==0?"+":"-")<<_i->chrn<<":"<<_i->ref_start_pos<<"-"<<_i->ref_end_pos<<"/"<<_i->qry_name<<":"<<_i->qry_start_pos<<"-"<<_i->qry_end_pos<<") "<<std::endl;
               }
            }
            if (_gp_pair_map.find(_fg_it->first)==_gp_pair_map.end()){
               _gp_pair_map[_fg_it->first] = oss_tostr.str();
            }else{
               std::cout<<"Error!!! one pair occur twice: "<<_fg_it->first<<std::endl;
               std::cout<<_gp_pair_map.find(_fg_it->first)->second;
               std::cout<<oss_tostr.str();
            }
            NumStrPair _snp;
            _snp._num_ = _d_dp_it->second.cand1.size();
            _snp._str_= _fg_it->first;
            //_snp._str_ = oss_tostr.str();
            _gp_rank.push_back(_snp);
 
         }
      }
   }
   std::map<std::string, std::vector<std::string> > _gf_time;
   std::map<std::string, std::vector<std::string> >::iterator _gf_t_i;
   for(std::map<std::string, std::string>::iterator _gpp_m_i=_gp_pair_map.begin(); _gpp_m_i!=_gp_pair_map.end(); _gpp_m_i++){
      _gf_t_i = _gf_time.find(_fg_cand[_gpp_m_i->first].g1);
      if (_gf_t_i==_gf_time.end()){
         _gf_time[_fg_cand[_gpp_m_i->first].g1] = std::vector<std::string>();
      }
      _gf_time[_fg_cand[_gpp_m_i->first].g1].push_back(_gpp_m_i->first);
      _gf_t_i = _gf_time.find(_fg_cand[_gpp_m_i->first].g2);
      if (_gf_t_i==_gf_time.end()){
         _gf_time[_fg_cand[_gpp_m_i->first].g2] = std::vector<std::string>();
      }
      _gf_time[_fg_cand[_gpp_m_i->first].g2].push_back(_gpp_m_i->first);
   }

   std::sort(_gp_rank.begin(), _gp_rank.end(), compare_support);
   for(std::vector<NumStrPair>::reverse_iterator _gp_it=_gp_rank.rbegin(); _gp_it!=_gp_rank.rend(); _gp_it++){
      if (_gf_time[_fg_cand[_gp_it->_str_].g2].size()>2 || _gf_time[_fg_cand[_gp_it->_str_].g1].size()>2) { continue; }
      else if ((_gf_time[_fg_cand[_gp_it->_str_].g2].size()>1 || _gf_time[_fg_cand[_gp_it->_str_].g1].size()>1) && _gp_it->_num_<10){ continue; }
      std::cout<<_gp_pair_map[_gp_it->_str_];
   }

   bam_destroy1(bam_one_alignment);
   bam_hdr_destroy(hdr);
   sam_close(in_bam);
   return 0; 
}



int usage(FILE * fp)
{
    fprintf (fp, "Usage: ./get_gfFrombam <input_bam> <input_gtf> <min-overlap-len> <bin_size> <min-map-len> [pseudogene:0(default)/1] [Secondary_alignment:0(default)] [min_sup_read:2(default)]\n");
    return 0;
}

int main (int argc, char * argv[])
{
   char * in_bam = NULL;
   char * in_gtf_file = NULL;
   int min_ovlp_len;
   bool has_error = false;
   int64_t _bin_size = 200;
   int64_t _min_map_len;
   int _used_pseudogene = 0;
   int _used_secondary_alignment = 0;
   int _min_sup_read = 2;

   if (argc<6){
      usage(stderr);
      return 1;
   }

   in_bam = argv[1];
   if (argc>7){
      fprintf(stdout, "%s %s %s %s %s %s\n", argv[1], argv[2], argv[3], argv[4],  argv[5], argv[6]);
   }else{
      fprintf(stdout, "%s %s %s %s %s\n", argv[1], argv[2], argv[3], argv[4], argv[5]);
   }
   in_gtf_file = argv[2];
   min_ovlp_len = atoi(argv[3]);
   _bin_size = atoi(argv[4]);
   _min_map_len = atoi(argv[5]);
   if (argc>6){
      _used_pseudogene = atoi(argv[6]);
   }
   if (argc>7){
      _used_secondary_alignment = atoi(argv[7]);
   }
   if (argc>8){
      _min_sup_read = atoi(argv[8]);
   }

   if (access(in_bam, F_OK)==-1){
      fprintf(stderr, "Bam file (%s) does not exist \n", in_bam);
      has_error = true;
   }

   if (access(in_gtf_file, F_OK)==-1){
      fprintf(stderr, "Bam file (%s) does not exist \n", in_gtf_file);
      has_error = true;
   }
   /*
   if (min_ovlp_len>100){
      fprintf(stderr, "Input max ovlp (%d) len between genes is larger than 100 \n", min_ovlp_len);
      has_error = true;
   }

   if (min_ovlp_len>50){
      fprintf(stderr, "Warning!! Input max ovlp (%d) len between genes is larger than 50 \n", min_ovlp_len);
   }
   */
   if (has_error){
      usage(stderr);
      return 1;
   }

   m_check_gene_fusion(in_bam, in_gtf_file, min_ovlp_len, _bin_size, _used_pseudogene, _min_map_len, _used_secondary_alignment, _min_sup_read);
}

