#ifndef _COM_FUN_H
#define _COM_FUN_H

#include "_com_structs_.h"

extern char m_cigar_str[]; // = "MIDNSHP=XB";

bool compare_support(const NumStrPair & snp1, const NumStrPair & snp2);

std::vector<std::string> m_split_string(const std::string & m_str, char m_delimiter); 

std::map<std::string, std::map<std::string, GenomicRegion> > get_gene_from_gtf(const char* in_gtf_file, const int _used_pseudogene);

#endif

