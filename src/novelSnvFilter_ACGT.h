#include <cstdio>
#include <getopt.h>
#include <cstdlib>
#include <cstring>


struct parameters {
  char* var_f;
  char* mapping_f;
  char* type;
  unsigned int unique;
  unsigned int skipPileup;
  char* chr;
};

struct parameters* interface(struct parameters* param,int argc, char *argv[]);
void delete_param(struct parameters* param);
void usage(void);

const char* program_name;

struct parameters* interface(struct parameters* param, int argc, char *argv[]){

  program_name = argv[0];
  int c;     // the next argument
  int help = 0;

  if (argc < 2) {
    usage();
    exit(0);
  }

  param = new struct parameters;
  param->var_f = new char;
  param->mapping_f = new char;
  param->type = new char;
  param->chr = new char;
 
  const struct option long_options[] ={
    {"var",1,0, 'v'},
    {"mapping",1,0,'m'},
    {"type",1,0,'t'},
    {"unique",0,0,'u'},
    {"skipPileup",0,0,'s'},
    {"chr",1,0,'c'},
    {"help",0,0,'h'},
    {0, 0, 0, 0}
  };

  while (1) {

    int option_index = 0;
    c = getopt_long_only (argc, argv,"husv:m:t:c:",long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch(c) {
    case 0:
      break;
    case 'v':
      param->var_f = optarg;
      break;
    case 'm':
      param->mapping_f = optarg;
      break;
    case 't':
      param->type = optarg;
      break;
    case 'u':
      param->unique = 1;
      break;
    case 's':
      param->skipPileup = 1;
      break;
    case 'c':
      param->chr = optarg;
      break;
    case 'h':
      help = 1;
      break;
    case '?':
      help = 1;
      break;
    default:
      help = 1;
      break;
    }
  }

  if(help) {
    usage();
    delete_param(param);
    exit(0);
  }

  return param;
}

void usage()
{
  fprintf(stdout, "\nnovelSnvFilter.cpp, Copyright (C) 2020 Sun Ruping <ruping@umn.edu>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "Usage: %s options [inputfile] \n\n", program_name);
  fprintf(stdout, "-h --help                print the help message\n");
  fprintf(stdout, "-v --var     <filename>  sorted vcf file contains only single nucleotide variants to be checked.\n");
  fprintf(stdout, "-m --mapping <filename>  mapping_file (coordinates' sorted bam-file or file of bam-file names).\n");
  fprintf(stdout, "-q --unique              only calculate for uniquely mapped reads.\n");
  fprintf(stdout, "-q --skipPileup          skip piled up reads.\n");
  fprintf(stdout, "-c --chr     <prefix>    set to prefix when the chromosome names in bam files starting with \'prefix\', e.g., chr, Chr or CHR.\n");
  fprintf(stdout, "-t --type    <p/s>       under development, do not set at this moment\n");
  fprintf(stdout, "\n");
}


void delete_param(struct parameters* param)
{
  delete(param->var_f);
  delete(param->mapping_f);
  delete(param->type);
  delete(param->chr);
  delete(param);
}
