#ifndef _COM_STRUCT_H
#define _COM_STRUCT_H

#include <vector>
#include <string>
#include <cstdint>
#include <map>


typedef struct {
   int64_t start_pos;
   int64_t end_pos;
   std::string chrn;
} GenomicRegion;


typedef struct {
   int64_t _num_;
   std::string _str_;
} NumStrPair;

typedef struct {
   int64_t qry_start_pos;
   int64_t qry_end_pos;
   int64_t ref_start_pos;
   int64_t ref_end_pos;
   int map_strand;
   std::string chrn;
   int64_t other_pos;
   std::string qry_name;
} GenomicMapRegion;


typedef struct {
  std::string g1;
  std::string g2;
  int64_t bin_pos_1;
  int64_t bin_pos_2;
  std::map<std::string, std::vector<GenomicMapRegion> > cand1;
  std::map<std::string, std::vector<GenomicMapRegion> > cand2;
} GeneFusionCand;



#endif 


