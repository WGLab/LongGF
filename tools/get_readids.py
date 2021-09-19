#!/usr/bin/env python


import os,sys

if __name__=='__main__':
   if len(sys.argv)<2:
      print("Usage: \n")
      print("       python {} {} {}".format(sys.argv[0], "log-file", "optional(interest-of-fusion-list)"))
      print("For example: ")
      print("       python {} {}".format(sys.argv[0], "log/longgf_output.log"));
      print("       python {} {} {}".format(sys.argv[0], "log/longgf_output.log", "BCAS4:BCAS3"))
      print("       python {} {} {}".format(sys.argv[0], "log/longgf_output.log", "BCAS4:BCAS3/MGAT5:IGLC7"))
      sys.exit(1)

   log_file = sys.argv[1];
   if len(sys.argv)>2:
      m_gf_list = sys.argv[2].split('/')
   else:
      m_gf_list = None;

   t_readids = {}
   t_gf = ''
   with open(log_file, 'r') as mr:
      is_interest = False;
      while True:
         line = mr.readline();
         if not line: break;
         
         lsp = line.strip().split();
         if len(lsp)==0: continue;

         if (lsp[0] in ['GF']): 
            if (not m_gf_list==None) and len(t_readids)>0:
               print("{}\n\t{}".format( t_gf, '\n\t'.join( list(sorted( t_readids.keys() )) ) ))
            if ( m_gf_list==None or lsp[1] in m_gf_list):
               is_interest = True; 
               t_gf = lsp[1]
               if not m_gf_list==None: t_readids = {}
            else: is_interest = False; t_readids = {}
         else:
            if is_interest and line[0] in [' ','\t']:
               t_readid = line.strip().split();
               if len(t_readid)>1: 
                  t_readid = t_readid[0].split('/')
                  if len(t_readid)>1:
                     t_readid = t_readid[1].split(':')
                  if len(t_readid)>1:
                     t_readids[ t_readid[0] ] = True;
      if len(t_readids)>0:
         if m_gf_list==None:
            print("{}".format( '\n'.join( list(sorted( t_readids.keys() )) ) ))   
         else:
            print("{}\n\t{}".format( t_gf, '\n\t'.join( list(sorted( t_readids.keys() )) ) ))


