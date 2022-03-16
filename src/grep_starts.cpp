/*****************************************************************************

  (c) 2020 - Sun Ruping
  ruping@umn.edu
  grep starts for sequencing reads

g++ grep_starts.cpp -I/home/ruping/ruping/tools/bamtools/include/bamtools/ -I/home/ruping/ruping/tools/zlib/current/include/ -I/home/ruping/ruping/tools/boost/current/include/ -L/home/ruping/ruping/tools/bamtools/lib64/ -L/home/ruping/ruping/tools/zlib/current/lib/ -L/home/ruping/ruping/tools/boost/current/lib/ -lbamtools -lz -Wl,-rpath,/home/ruping/ruping/tools/bamtools/lib64/:/home/ruping/ruping/tools/boost/current/lib/ -lboost_regex -pthread -static -o grep_starts

g++ grep_starts.cpp 
-I/home/regularhand/tools/bamtools/include/ -I/home/regularhand/tools/zlib/current/include/ -I/home/regularhand/tools/boost/current/include/ 
-L/home/regularhand/tools/bamtools/lib/ -L/home/regularhand/tools/zlib/current/lib/ -L/home/regularhand/tools/boost/current/lib/ 
-lbamtools -lz -Wl,-rpath,/home/regularhand/tools/bamtools/lib/:/home/regularhand/tools/boost/current/lib/ -lboost_regex -o grep_starts

******************************************************************************/

#include <api/BamReader.h>
#include <api/BamMultiReader.h>
using namespace BamTools;

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <deque>
#include <set>
#include <string>
#include <cstring>
#include <sstream>
#include "grep_starts.h"
using namespace std;

struct region {  // a bed file containing gene annotations
  string chr;
  string chro;
  unsigned int start;
  unsigned int end;
  //string name;
  //float score;
  //string strand;
  //unsigned int thickStart;
  //unsigned int thickEnd;
  //string itemRgb;
  //unsigned int blockCount;
  //string blockSizes;
  //string blockStarts;
  // results storing here
  unsigned int tags;  // storing position and count of tags relative to a gene
  unsigned int starts;
};


unsigned int read_length = 0;

inline void ParseCigar(const vector<CigarOp> &cigar, vector<int> &blockStarts, vector<int> &blockEnds, unsigned int &alignmentEnd);
inline void splitstring(const string &str, vector<string> &elements, const string &delimiter);
inline void eatline(const string &str, deque <struct region> &region_ref, bool &withChr);
inline string int2str(unsigned int &i);
inline string float2str(float &f);
inline void gene_processing(struct region &gene);

int main ( int argc, char *argv[] ) { 

  struct parameters *param = 0;
  param = interface(param, argc, argv);

  //region file input (the region file should be sorted as the same way as the bam file)
  ifstream region_f;
  region_f.open(param->region_f, ios_base::in);  // the region file is opened


  //bam input and generate index if not yet 
  //-------------------------------------------------------------------------------------------------------+
  // BAM input (file or filenames?)                                                                        |
  //-------------------------------------------------------------------------------------------------------+
  char *fof = param->mapping_f;
  FILE *IN=NULL;
  char linefof[5000];
  int filecount=0;
  vector <string> fnames;

  if (strchr(fof,' ')!=NULL) {
    char *ptr;
    ptr=strtok(fof," ");
    while (ptr!=NULL) {
      fnames.push_back(ptr);
      filecount++;
      ptr=strtok(NULL," ");
    }
  } else {
    IN=fopen(fof,"rt");
    if (IN!=NULL) {
      long linecount=0;
      while (fgets(linefof,5000-1,IN)!=NULL) {
        linecount++;
        if (linefof[0]!='#' && linefof[0]!='\n') {
          char *ptr=strchr(linefof,'\n');
          if (ptr!=NULL && ptr[0]=='\n') {
            ptr[0]='\0';
          }
          FILE *dummy=NULL;
          dummy=fopen(linefof,"rt");
          if (dummy!=NULL) {     // seems to be a file of filenames...
            fclose(dummy);
            fnames.push_back(linefof);
            filecount++;
          } else if (filecount==0 || linecount>=1000-1) {  // seems to be a single file
            fnames.push_back(fof);
            filecount++;
            break;
          }
        }
      }
      fclose(IN);
    }
  }  //file or file name decided and stored in vector "fnames"

  cerr << "the input mapping files are:" << endl;
  vector <string>::iterator fit = fnames.begin();
  for(; fit != fnames.end(); fit++) {
    cerr << *fit << endl;
  }

  //-------------------------------------------------------------------------------------------------------+
  // end of file or filenames                                                                              |
  //-------------------------------------------------------------------------------------------------------+

  // open the BAM file(s)
  BamMultiReader reader;
  reader.Open(fnames);

  // get header & reference information
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();

  if ( ! reader.LocateIndexes() )     // opens any existing index files that match our BAM files
    reader.CreateIndexes();         // creates index files for BAM files that still lack one


  //should decide which chromosome
  string line;
  string old_chr = "SRP";
  string type = param->type;
  bool startwithChr = false;
  if (param->chr == 1){
    startwithChr = true;
  }

  //regions for the input of region file
  deque <struct region> regions;

  getline(region_f, line); //get the first line
  eatline(line,regions, startwithChr);
  
  deque <struct region>::iterator it = regions.begin();

  while ( it->chr != old_chr ) {

    old_chr = it->chr;  // set the current chr as old chr

    int chr_id  = reader.GetReferenceID(it->chr);

    if ( chr_id == -1 ) {  //reference not found

      for (; it != regions.end() && it->chr == old_chr; ) {
        gene_processing(*it);           // print the old region info
        it = regions.erase(it);         // erase the current region
      }
  
      while ( regions.empty() ) {    
        getline(region_f, line);
        if ( region_f.eof() ){
          cerr << "finished: end of region file, zone 0" << endl;
          break;
        }
        eatline(line, regions, startwithChr);
        it = regions.begin();
        if (it->chr == old_chr){  
          gene_processing(*it);      
          regions.clear();
          continue;
        }
      }
      continue;
    }

    int chr_len = refs.at(chr_id).RefLength;

    if ( !reader.SetRegion(chr_id, 1, chr_id, chr_len) ) // here set region
      {
        cerr << "bamtools count ERROR: Jump region failed " << it->chr << endl;
        reader.Close();
        exit(1);
      }

    BamAlignment bam;
    while (reader.GetNextAlignment(bam)) {

      if ( bam.IsMapped() == false ) continue;              // skip unaligned reads
      if ( bam.IsDuplicate() == true ) continue;            // skip PCR duplicates

      unsigned int unique = 0;
      //if ( bam.HasTag("NH") ) {
      // bam.GetTag("NH", unique);                   // uniqueness
      //} else if (bam.HasTag("XT")) {
      //  string xt;
      //  bam.GetTag("XT", xt);                       // bwa aligner
      //  xt = xt.substr(0,1);
      //  if (xt != "R") {
      //    unique = 1;
      //  }
      //} else {
        if (bam.MapQuality > 10 || bam.MapQuality == 0) {                   // bowtie2
          unique = 1;
        }
        //}

      if (param->unique == 1) {
        if (unique != 1) {                         // skipe uniquelly mapped reads              
          continue;
        }
      }


      if (read_length == 0){
        read_length = bam.Length;
      }

      string chrom = refs.at(bam.RefID).RefName;
      string strand = "+";
      if (bam.IsReverseStrand()) strand = "-";

      unsigned int alignmentStart =  bam.Position+1;
      //unsigned int mateStart = bam.MatePosition+1;
      unsigned int alignmentEnd = bam.GetEndPosition(false, true);
      unsigned int cigarEnd;
      vector <int> blockLengths;
      vector <int> blockStarts;
      blockStarts.push_back(0);
      ParseCigar(bam.CigarData, blockStarts, blockLengths, cigarEnd);


      deque <struct region>::iterator iter = regions.begin();

      if ( iter->start > alignmentEnd ) continue;  // skip reads not overlapping with the first region

      while ( iter->chr == old_chr && iter->start <= alignmentEnd && iter != regions.end() ) {

        if (iter->end < alignmentStart) {            // the region end is beyond the alignmentStart

          gene_processing(*iter);                    // processing
          iter = regions.erase(iter);                // this region should be removed
          if ( regions.empty() ) { 
            getline(region_f, line);                        // get a line of region file
            if ( ! region_f.eof() ) {
              eatline(line, regions, startwithChr);                         // eat a line and put it into the duque
              iter = regions.begin();
            }
            else {  // it's reaching the end of the region file
              cerr << "finished: end of region file, zone 3" << endl;
              break;
            }
          }
          continue;
        }

        if (iter->end >= alignmentStart && iter->start <= alignmentEnd) {  //overlapping, should take action

           iter->tags += 1;

           if (alignmentStart >= iter->start && alignmentStart <= iter->end) {
             iter->starts += 1;
           }

        }  // overlapping take action!

        if ( (iter+1) != regions.end() )
          iter++;                                           // if this region is not the last element in the deque
        else {                                              // the last element
          getline(region_f, line);                          // get a line of region file
          if ( ! region_f.eof() ){
            eatline(line, regions, startwithChr);                         // eat a line and put it into the duque
            iter = regions.end();
            iter--;
          }
          else {  //it's reaching the end of the region file
            cerr << "finished: end of region file, zone 4" << endl;
            break;
          }
        }

      } //while

    }  // read a bam

 
    //somehow to loop back
    it = regions.begin();                   //reset to begin
    for (; it != regions.end() && it->chr == old_chr; ) {
      gene_processing(*it);              // print the old region info
      it = regions.erase(it);             // erase the current region
    }
  
    while ( regions.empty() ) {    

      getline(region_f, line);
      if ( region_f.eof() ){
        cerr << "finished: end of region file, zone 5" << endl;
        exit(0);
      }
      eatline(line, regions, startwithChr);
      it = regions.begin();
      if (it->chr == old_chr){
        gene_processing(*it);
        regions.clear();
        continue;
      }
    }

  } // region chr != old chr
      
  regions.clear();
  reader.Close();
  region_f.close();
  return 0;

} //main


inline string int2str(unsigned int &i){
  string s;
  stringstream ss(s);
  ss << i;
  return ss.str();
}


inline string float2str(float &f){
  string s;
  stringstream ss(s);
  ss << f;
  return ss.str();
}


inline void splitstring(const string &str, vector<string> &elements, const string &delimiter) {
  string::size_type lastPos = str.find_first_not_of(delimiter, 0);
  string::size_type pos     = str.find_first_of(delimiter, lastPos);

  while (string::npos != pos || string::npos != lastPos) {
    elements.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiter, pos);
    pos = str.find_first_of(delimiter, lastPos);
  }
}


inline void eatline(const string &str, deque <struct region> &region_ref, bool &withChr) {
  
   vector <string> line_content;
   //split line and then put it into a deque

   splitstring(str, line_content, "\t");
   vector <string>::iterator iter = line_content.begin();
   unsigned int i;

   struct region tmp;
   tmp.tags = 0;
   tmp.starts = 0;
 
   for(i = 1; iter != line_content.end(); iter++, i++){
     switch (i) {
     case 1:  // chr
       tmp.chr = *iter;
       tmp.chro = *iter;
       if (withChr == true) {
         if ((tmp.chr).substr(0,1) != "c" && (tmp.chr).length() < 3) {   //mostlikely not starting with chr
           tmp.chr = "chr"+tmp.chr;
         }
         if(tmp.chr == "chrMT"){
           tmp.chr = "chrM";
         }
       } else { //'chr' is not required
         if ((tmp.chr).substr(0,1) == "c" && (tmp.chr).length() > 3) {   //mostlikely starting with chr
           tmp.chr = (tmp.chr).substr(3);
         }
       }
       continue;
     case 2:  // start
       tmp.start = atoi((*iter).c_str()) + 1;
       continue;
     case 3:  // end
       tmp.end = atoi((*iter).c_str());
       continue;
       //case 4:  // name
       //tmp.name = *iter;
       //continue;
       //case 5:  // score
       //tmp.score = atof((*iter).c_str());
       //continue;
       //case 6:  // strand
       //tmp.strand = *iter;
       //continue;
       //case 7:  // thickStart
       //tmp.thickStart = atoi((*iter).c_str());
       //continue;
       //case 8:  // thickEnd
       //tmp.thickEnd = atoi((*iter).c_str());
       //continue;
       //case 9: // itemRgb
       //tmp.itemRgb = *iter;
       //continue;
       //case 10: // blockCount
       //tmp.blockCount = atoi((*iter).c_str());
       //continue;
       //case 11: // blockSizes
       //tmp.blockSizes = *iter;
       //continue;
       //case 12: // blockStarts
       //tmp.blockStarts = *iter;
       //continue;
     default:
       break;
     }
   }

   region_ref.push_back(tmp);
}


inline void ParseCigar(const vector<CigarOp> &cigar, vector<int> &blockStarts, vector<int> &blockLengths, unsigned int &alignmentEnd) {

  int currPosition = 0;
  int blockLength  = 0;

  //  Rip through the CIGAR ops and figure out if there is more
  //  than one block for this alignment
  vector<CigarOp>::const_iterator cigItr = cigar.begin();
  vector<CigarOp>::const_iterator cigEnd = cigar.end();
  for (; cigItr != cigEnd; ++cigItr) {
    switch (cigItr->Type) {
    case ('M') :                           // matching
      blockLength  += cigItr->Length;
      currPosition += cigItr->Length;
      break;
    case ('I') : break;                    // insertion
    case ('S') : break;                    // soft-clipping
    case ('D') : break;                    // deletion
      blockLength  += cigItr->Length;
      currPosition += cigItr->Length;
      break;
    case ('P') : break;                    // padding
    case ('N') :                           // skipped region
      blockStarts.push_back(currPosition + cigItr->Length);
      blockLengths.push_back(blockLength);
      currPosition += cigItr->Length;
      blockLength = 0;                     // a new block
      break;
    case ('H') : break;                    // for 'H' - do nothing, move to next op
    default    :
      printf("ERROR: Invalid Cigar op type\n");   // shouldn't get here
      exit(1);
    }
  }
  // add the kast block and set the
  // alignment end (i.e., relative to the start)
  blockLengths.push_back(blockLength);
  alignmentEnd = currPosition;
}


inline void gene_processing(struct region &gene) {

  cout << gene.chro << "\t" << gene.start << "\t" << gene.end << "\t" << gene.tags << "\t" << gene.starts << endl;

}


