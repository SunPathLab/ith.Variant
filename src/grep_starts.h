#include <cstdio>
#include <getopt.h>
#include <cstdlib>
#include <cstring>


struct parameters {
  char* region_f;
  char* mapping_f;
  char* type;
  unsigned int unique;
  unsigned int chr;
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
  param->region_f = new char;
  param->mapping_f = new char;
  param->type = new char;
 
  const struct option long_options[] ={
    {"region",1,0, 'r'},
    {"mapping",1,0,'m'},
    {"type",1,0,'t'},
    {"unique",0,0,'u'},
    {"chr",0,0,'c'},
    {"help",0,0,'h'},
    {0, 0, 0, 0}
  };

  while (1) {

    int option_index = 0;
    c = getopt_long_only (argc, argv,"hur:m:t:c",long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch(c) {
    case 0:
      break;
    case 'r':
      param->region_f = optarg;
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
    case 'c':
      param->chr = 1;
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
  fprintf(stdout, "\ngrep_starts, Copyright (C) 2016 Sun Ruping <ruping@stanford.edu>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "Usage: %s options [inputfile] \n\n", program_name);
  fprintf(stdout, "-h --help    print the help message\n");
  fprintf(stdout, "-r --region  <filename>  UCSC gene annotation file in 12 column bed format (should be sorted according to chromosomes and coordinates)\n");
  fprintf(stdout, "-m --mapping <filename>  mapping_file (RNA-seq bam file, chromosomes and coordinates sorted also)\n");
  fprintf(stdout, "-q --unique              only calculate for uniquely mapped reads (set this when the bam files contain multi-mapping reads).\n");
  fprintf(stdout, "-c --chr                 set when the chromosome names in bam files starting with \'chr\'.\n");
  fprintf(stdout, "-t --type    <p/s>       under development\n");
  fprintf(stdout, "\n");
}


void delete_param(struct parameters* param)
{
  delete(param->region_f);
  delete(param->mapping_f);
  delete(param->type);
  delete(param);
}
