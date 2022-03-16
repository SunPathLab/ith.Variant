/*****************************************************************************

  (c) 2020 - Sun Ruping
  ruping@umn.edu

g++ novelSnvFilter_ACGT.cpp -I/home/ruping/ruping/tools/bamtools/include/bamtools/ -I/home/ruping/ruping/tools/zlib/current/include/ -I/home/ruping/ruping/tools/boost/current/include/ -L/home/ruping/ruping/tools/bamtools/lib64/ -L/home/ruping/ruping/tools/zlib/current/lib/ -L/home/ruping/ruping/tools/boost/current/lib/ -lbamtools -lz -Wl,-rpath,/home/ruping/ruping/tools/bamtools/lib64/:/home/ruping/ruping/tools/boost/current/lib/ -lboost_regex -pthread -static -o novelSnvFilter_ACGT

g++ novelSnvFilter_ACGT.cpp
-I/home/regularhand/tools/bamtools/include/ -I/home/regularhand/tools/zlib/current/include/ -I/home/regularhand/tools/boost/current/include/ 
-L/home/regularhand/tools/bamtools/lib/ -L/home/regularhand/tools/zlib/current/lib/ -L/home/regularhand/tools/boost/current/lib/ 
-lbamtools -lz -Wl,-rpath,/home/regularhand/tools/bamtools/lib/:/home/regularhand/tools/boost/current/lib/ -lboost_regex -o novelSnvFilter_ACGT
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
#include <iomanip>
#include "boost/regex.hpp"
#include "novelSnvFilter_ACGT.h"
using namespace std;
using namespace boost;


struct var {  // a bed file containing gene annotations
  string chr;
  string chro;  //original chr name from the input list
  string snpID;
  string ref;
  string alt;
  unsigned int start;
  unsigned int end;
  // results storing here
  unsigned int countAlt;
  unsigned int countAll;
  unsigned int countA;
  unsigned int countAn;
  unsigned int countC;
  unsigned int countCn;
  unsigned int countG;
  unsigned int countGn;
  unsigned int countT;
  unsigned int countTn;
  unsigned int countMappingGood;
  unsigned int countMappingBad;
  unsigned int countPositive;
  unsigned int countNegative;
  unsigned int inends;
  unsigned int countJump;
  unsigned int F1R2_alt;
  unsigned int F1R2_all;
  unsigned int F2R1_alt;
  unsigned int F2R1_all;
  unsigned int readlen;
  vector <unsigned int> lenVarReads;
  vector <unsigned int> surrounding;
  vector <unsigned int> surroundingIndels;
  map <unsigned int, unsigned int> conMis;
  string qualities;
};


//unsigned int read_length = 0;

inline void ParseCigar(const vector<CigarOp> &cigar, vector<int> &blockStarts, vector<int> &blockEnds, unsigned int &alignmentEnd, map<unsigned int, unsigned int> &insertions, unsigned int &softClip);
inline void splitstring(const string &str, vector<string> &elements, const string &delimiter);
inline bool eatline(const string &str, deque <struct var> &var_ref, string &withChr);
inline string int2str(unsigned int &i);
inline string float2str(float &f);
inline void var_processing(struct var &variant);
inline float CalcMedian (vector<unsigned int> &scores);

int main ( int argc, char *argv[] ) {

  struct parameters *param = 0;
  param = interface(param, argc, argv);

  //region file input (the region file should be sorted as the same way as the bam file)
  ifstream var_f;
  var_f.open(param->var_f, ios_base::in);  // the region file is opened


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
  string startwithChr = param->chr;
  if (startwithChr == "") {
    startwithChr = "none";
  }
  cerr << "chr prefix is: " << startwithChr << endl;

  //regions for the input of region file
  deque <struct var> variants;

  bool isComment;
  getline(var_f, line); //get the first line
  isComment = eatline(line, variants, startwithChr);
  while (isComment == true){
    //print comment
    //cout << line << endl;
    getline(var_f, line); //get a new line
    isComment = eatline(line, variants, startwithChr); //eat a new line
  }
  
  deque <struct var>::iterator it = variants.begin();

  while ( it->chr != old_chr ) {

    old_chr = it->chr;  // set the current chr as old chr

    int chr_id  = reader.GetReferenceID(it->chr);

    if ( chr_id == -1 ) {  //reference not found

      for (; it != variants.end() && it->chr == old_chr; ) {
        var_processing(*it);           // print the old region info
        it = variants.erase(it);         // erase the current region
      }
  
      while ( variants.empty() ) {    
        getline(var_f, line);
        if ( var_f.eof() ){
          cerr << "finished: end of region file, zone 0" << endl;
          break;
        }
        eatline(line, variants, startwithChr);
        it = variants.begin();
        if (it->chr == old_chr){  
          var_processing(*it);      
          variants.clear();
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

    ////pile-up pos stats   
    //map <string, unsigned int> PileUp;
    //unsigned int old_PileUp_start = 0;

    BamAlignment bam;
    while (reader.GetNextAlignment(bam)) {

      if ( bam.IsMapped() == false ) continue;      // skip unaligned reads
      if ( bam.IsDuplicate() == true && param->skipPileup == 1) continue;            // skip PCR duplicates

      unsigned int unique = 0;
      if ( bam.HasTag("NH") ) {
        bam.GetTag("NH", unique);                   // uniqueness
      } else {
        if (bam.MapQuality > 10) {                  // other aligner
          unique = 1;
        }
      }

      if (param->unique == 1) {
        if (unique != 1) {                         // skipe uniquelly mapped reads
          continue;
        }
      }


      //if (bam.Length > read_length) {              // get the read length
      //  read_length = bam.Length;
      //}

      string chrom = refs.at(bam.RefID).RefName;
      string strand = "+";
      if (bam.IsReverseStrand()) strand = "-";
      string FxRx = "FxRx";
      if (bam.IsProperPair()) {
        if (bam.IsFirstMate()) {   //first mate
          if (strand == "+") {
            FxRx = "F1R2";
          } else {
            FxRx = "F2R1";
          }
        } else {                   //second mate
          if (strand == "+") {
            FxRx = "F2R1";
          } else {
            FxRx = "F1R2";
          }
        }
      }
      unsigned int mappingQuality = bam.MapQuality;

      unsigned int alignmentStart =  bam.Position+1;
      unsigned int alignmentEnd = bam.GetEndPosition();

      unsigned int cigarEnd;
      vector <int> blockLengths;
      vector <int> blockStarts;
      map<unsigned int, unsigned int> insertions;       // for insertions 
      unsigned int softClip = 0;                        // for soft clipping
      blockStarts.push_back(0);
      ParseCigar(bam.CigarData, blockStarts, blockLengths, cigarEnd, insertions, softClip);


      //// do pileup check for duplicates
      //string alignSum = int2str(alignmentStart) + "\t" + bam.QueryBases;
      //if ( alignmentStart != old_PileUp_start ) {
      //  PileUp.clear();           //clear PileUp set                                                                                                                                                                                             
      //  PileUp.insert( pair <string, unsigned int> (alignSum, 1) );  //insert the new read
      //}  else if ( alignmentStart == old_PileUp_start ) { // same starts
      //  if ( PileUp.count(alignSum) > 0 ) {  // PileUp                                                   
      //    PileUp[alignSum]++;
      //    if ( bam.IsDuplicate() == true ) {            // skip PCR duplicates
      //      continue;
      //    } //PCR duplicates                       
      //  } else {
      //    PileUp.insert( pair <string, unsigned int> (alignSum, 1) );
      //  }
      //} //same starts                                                      
      //old_PileUp_start = alignmentStart;
      ////pile up check


      deque <struct var>::iterator iter = variants.begin();

      if ( iter->start > alignmentEnd ) continue;          // skip reads not overlapping with the first region

      while ( iter->chr == old_chr && iter->start <= alignmentEnd && iter != variants.end() ) {

        if (iter->end < alignmentStart) {                  // the region end is beyond the alignmentStart

          var_processing(*iter);                           // processing
          iter = variants.erase(iter);                     // this region should be removed
          if ( variants.empty() ) { 
            getline(var_f, line);                          // get a line of region file
            if ( ! var_f.eof() ) {
              eatline(line, variants, startwithChr);                     // eat a line and put it into the duque
              iter = variants.begin();
            }
            else {  // it's reaching the end of the region file
              cerr << "finished: end of region file, zone 3" << endl;
              break;
            }
          }
          continue;
        }

        if ( iter->end >= alignmentStart && iter->start <= alignmentEnd ) {  //overlapping, should take action

          if (bam.Length > iter->readlen) {                 // should we re-define the read length?
            iter->readlen = bam.Length;
          }
          
          unsigned int mismatches = 0;                      // how many mismatches (including indels) does this read have?
          unsigned int indels = 0;                          // how many indels does this read have?
          bool varInRead = false;                           // is the var in the read?
          bool posInRead = false;
          vector <int>::iterator bliter = blockLengths.begin();
          vector <int>::iterator bSiter = blockStarts.begin();
          while (bliter != blockLengths.end() && bSiter != blockStarts.end()) {
            unsigned int blockstart = *bSiter + alignmentStart;
            unsigned int blockend = *bliter + blockstart;
            if (iter->start >= blockstart && iter->end <= blockend) {
               posInRead = true;
               break;
            } //overlap
            bliter++;
            bSiter++;
          }

          if (posInRead == true) {    //need to get strand information for all reads !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             iter->countAll += 1;
             if (strand == "+") {
               iter->countPositive += 1;
             } else {
               iter->countNegative += 1;
             }
             if (FxRx == "F1R2") {
               iter->F1R2_all += 1;
             } else if (FxRx == "F2R1") {
               iter->F2R1_all += 1;
             }
          }
          else
            iter->countJump += 1;

          //processing MD string, calculate mismatch coordinates and compare with the variants 
          string MD;
          bam.GetTag("MD", MD);
          //cerr << MD << endl;

          unsigned int cuPos = alignmentStart;
          unsigned int cuPosRead = softClip + 1;

          //cerr << iter->start << "\t" << bam.Name << "\t" << cuPosRead << endl;      //deBUG!!!!!!!!

          regex rgx( "([0-9]+)([ACGTacgt]|\\^[ACGTacgt]+)" );
          int subs[] = {1,2};
          sregex_token_iterator rit ( MD.begin(), MD.end(), rgx, subs );
          sregex_token_iterator rend;

          map<unsigned int, unsigned int>::iterator inserit_index = insertions.begin();
          while ( inserit_index != insertions.end() ) {    // check insertions
            mismatches += 1;                               // count as mismatches
            indels += 1;                                   // count for indels
            inserit_index++;
          }
          inserit_index = insertions.begin();              //reset it for the begin of insertions


          while ( rit != rend ) {

            unsigned int incre = atoi((*rit).str().c_str());                  //number 1
            cuPos += incre;                                                   //number 1
            cuPosRead += incre; 

            if (blockStarts.size() > 1) {                 //judge which block the mutation locate
              vector <int>::iterator bliter2 = blockLengths.begin();
              vector <int>::iterator bSiter2 = blockStarts.begin();
              unsigned int culength = 0;
              while (bliter2 != blockLengths.end() && bSiter2 != blockStarts.end()) {
                if (cuPosRead <= (culength + *bliter2)) {
                  cuPos += (*bSiter2 - culength);
                  break;
                }
                culength += *bliter2;
                bliter2++;
                bSiter2++;
              }               
            } //multi blocks especially useful for RNA-seq junction reads

            map<unsigned int, unsigned int>::iterator inserit = inserit_index;
            while ( inserit != insertions.end() ) {
              if ( inserit->first < cuPosRead ) {
                cuPosRead += inserit->second;
                inserit++;
                inserit_index = inserit;
              } else {
                inserit_index = inserit;
                break;
              }
            }

            ++rit;                                            //round 1 addition

            if (((*rit).str())[0] == '^') {                   //variant 2
              incre = (*rit).length() - 1;                    //variant 2
              cuPos += incre;                                 //variant 2
              mismatches += 1;
              indels += 1;
            } else if ((*rit).length() == 1) {                // single base nucleotide change

              //check whether it is "N" or not
              string baseInReadPre = (bam.QueryBases).substr( cuPosRead-1, 1 );
              if (baseInReadPre != "N") {
                 mismatches += 1;
                 map<unsigned int, unsigned int>::iterator cmi = (iter->conMis).find(cuPos);
                 if ( cmi == (iter->conMis).end() ) {                                            // not found need to record mismatch in a map
                   (iter->conMis).insert( pair <unsigned int, unsigned int> (cuPos, 1) );        // not found need to record mismatch in a map
                 } else {
                   cmi->second += 1;                                                             // found increase it
                 }
              }

              if ( cuPos == iter->start ) { // it is right here with some variant base!!!

                varInRead = true;
                
                if ((alignmentEnd - cuPos) <= 10 || (cuPos - alignmentStart) <= 10) {        // inends
                  iter->inends += 1;
                }

                if ( mappingQuality >= 30 ) {         //good mapping qual
                    iter->countMappingGood += 1;
                } else if (mappingQuality <= 29) {    // bad mapping qual
                    iter->countMappingBad += 1;
                }

                //cout << cuPosRead << endl;             // deBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                string baseInRead = (bam.QueryBases).substr( cuPosRead-1, 1 );
                if (baseInRead == iter->alt) {           // it is exactly the same alt base
                  iter->qualities += (bam.Qualities).substr( cuPosRead-1, 1 );   //base quality
                  if (FxRx == "F1R2") {
                    iter->F1R2_alt += 1;
                  } else if (FxRx == "F2R1") {
                    iter->F2R1_alt += 1;
                  }
                }
                
                iter->countAlt += 1;
                if (strand == "+") {                  //positive strand
                  if (baseInRead == "A") {
                    iter->countA += 1;
                  } else if (baseInRead == "C") {
                    iter->countC += 1;
                  } else if (baseInRead == "G") {
                    iter->countG += 1;
                  } else if (baseInRead == "T") {
                    iter->countT += 1;
                  }
                } else {                              //negative strand
                  if (baseInRead == "A") {
                    iter->countAn += 1;
                  } else if (baseInRead == "C") {
                    iter->countCn += 1;
                  } else if (baseInRead == "G") {
                    iter->countGn += 1;
                  } else if (baseInRead == "T") {
                    iter->countTn += 1;
                  }
                }
              }
              cuPos += 1;
              cuPosRead += 1;
            } else {
              cerr << "wired thing happened in the MD string of " << bam.Name << endl;
              exit(1);
            }

            ++rit;                                            //round 2 addition

          } //loop for all MD characters

          if (varInRead == true) {
            (iter->surrounding).push_back(mismatches);
            (iter->surroundingIndels).push_back(indels);
            (iter->lenVarReads).push_back(bam.Length);
          }

        }  // overlapping take action!

        if ( (iter+1) != variants.end() )
          iter++;                                           // if this region is not the last element in the deque
        else {                                              // the last element
          getline(var_f, line);                             // get a line of region file
          if ( ! var_f.eof() ) {
            eatline(line, variants, startwithChr);                        // eat a line and put it into the duque
            iter = variants.end();
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
    it = variants.begin();                   //reset to begin
    for (; it != variants.end() && it->chr == old_chr; ) {
      var_processing(*it);              // print the old region info
      it = variants.erase(it);             // erase the current region
    }
  
    while ( variants.empty() ) {    

      getline(var_f, line);
      if ( var_f.eof() ){
        cerr << "finished: end of region file, zone 5" << endl;
        exit(0);
      }
      eatline(line, variants, startwithChr);
      it = variants.begin();
      if (it->chr == old_chr){
        var_processing(*it);
        variants.clear();
        continue;
      }
    }

  } // region chr != old chr
      
  variants.clear();
  reader.Close();
  var_f.close();
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


inline bool eatline(const string &str, deque <struct var> &var_ref, string &withChr) {
  
  bool isComment = false;
  if (str[0] == '#' || str[0] == '@') {
    isComment = true;
    return isComment;
  }

  vector <string> line_content;
  //split line and then put it into a deque

  splitstring(str, line_content, "\t");
  vector <string>::iterator iter = line_content.begin();
  unsigned int i;

  struct var tmp;
  tmp.countAlt = 0;
  tmp.countAll = 0;
  tmp.countA = 0;
  tmp.countAn = 0;
  tmp.countC = 0;
  tmp.countCn = 0;
  tmp.countG = 0;
  tmp.countGn = 0;
  tmp.countT = 0;
  tmp.countTn = 0;
  tmp.countMappingGood = 0;
  tmp.countMappingBad = 0;
  tmp.countPositive = 0;
  tmp.countNegative = 0;
  tmp.inends = 0;
  tmp.countJump = 0;
  tmp.F1R2_alt = 0;
  tmp.F1R2_all = 0;
  tmp.F2R1_alt = 0;
  tmp.F2R1_all = 0;
  tmp.readlen = 0;
  
  for(i = 1; iter != line_content.end(); iter++, i++) {
    switch (i) {
    case 1:  // chr
      tmp.chr = *iter;
      tmp.chro = *iter;
      if (withChr != "none") {
        if ((tmp.chr).substr(0,1) != "c" && (tmp.chr).substr(0,1) != "C" && (tmp.chr).length() < 3) {   //mostlikely not starting with chr
          tmp.chr = withChr + tmp.chr;
        }
        //if(tmp.chr == "chrMT") {
        //  tmp.chr = "chrM";
        //}
      } else { //'chr' is not required
        if ( ( (tmp.chr).substr(0,1) == "c" || (tmp.chr).substr(0,1) == "C" ) && (tmp.chr).length() > 3) {   //mostlikely starting with chr
          tmp.chr = (tmp.chr).substr(3);
        }
      }
      continue;
    case 2:  // pos
      tmp.start = atoi((*iter).c_str());
      tmp.end = atoi((*iter).c_str());
      continue;
    case 3:  // id
      tmp.snpID = *iter;
      continue;
    case 4:  // ref
      tmp.ref = *iter;
      continue;
    case 5:  // alt
      tmp.alt = *iter;
      continue;
    default:
      break;
    }
  }

  var_ref.push_back(tmp);
  return isComment;

}


inline void ParseCigar(const vector<CigarOp> &cigar, vector<int> &blockStarts, vector<int> &blockLengths, unsigned int &alignmentEnd, map<unsigned int, unsigned int> &insertions, unsigned int &softClip) {

  int currPosition = 0;
  int blockLength  = 0;

  unsigned int insertSize;
  unsigned int insertPos;

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
    case ('I') :                           // insertion
      insertSize = cigItr->Length;
      insertPos = currPosition + 1;
      insertions.insert( pair <unsigned int, unsigned int> (insertPos, insertSize) );
      break;
    case ('S') :                           // soft-clipping
      if (currPosition == 0){ //only take action for the beginning clipping
        softClip = cigItr->Length;
      }
      break;
    case ('D') :                           // deletion
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


inline void var_processing(struct var &variant) {

  unsigned int ssum = 0;
  vector <unsigned int>::iterator sit = (variant.surrounding).begin();
  for(; sit != (variant.surrounding).end(); sit++) {
    ssum += *sit;
  }

  float meanMis;
  float medianMis;
  unsigned int surrSize = variant.surrounding.size();
  if (surrSize == 0) {
    meanMis = 0.0;
    medianMis = 0.0;
  } else {
    meanMis = ((float)ssum)/((float)surrSize);
    medianMis = CalcMedian(variant.surrounding);
  }

  unsigned int indelsum = 0;
  vector <unsigned int>::iterator indit = (variant.surroundingIndels).begin();
  for(; indit != (variant.surroundingIndels).end(); indit++) {
      indelsum += *indit;
  }

  float meanIndel;
  float medianIndel;
  unsigned int surrIndelSize = variant.surroundingIndels.size();
  if (surrIndelSize == 0) {
    meanIndel = 0.0;
    medianIndel = 0.0;
  } else {
    meanIndel = ((float)indelsum)/((float)surrIndelSize);
    medianIndel = CalcMedian(variant.surroundingIndels);
  }

  //unsigned int lenSum;
  //vector <unsigned int>::iterator lenit = (variant.lenVarReads).begin();
  //for(; lenit != (variant.lenVarReads).end(); lenit++) {
  //  lenSum += *lenit;
  //}
  
  //float meanVRLength;        // get Length of reads with variants
  float medianVRLength;      // get Length of reads with variants
  unsigned int lenVRSize = variant.lenVarReads.size();
  if (surrIndelSize == 0) {
    //meanVRLength = 0.0;
    medianVRLength = 0.0;
  } else {
    //meanVRLength = ((float)lenSum)/((float)lenVRSize);
    medianVRLength = CalcMedian(variant.lenVarReads);
  }
  

  float fracBadMappingQual = 0;
  if ((variant.countMappingGood + variant.countMappingBad) > 0) {
    fracBadMappingQual = ((float)(variant.countMappingBad))/((float)(variant.countMappingGood + variant.countMappingBad));
  }

  // get local error rate estimate
  float totalBases = (float)variant.countAll * (float)variant.readlen;
  map<unsigned int, unsigned int>::iterator cmi = (variant.conMis).begin();
  unsigned int numncMis = 0;
  for (; cmi != (variant.conMis).end(); cmi++) {
    if (cmi->second == 1) {
      numncMis += 1;
    }
  }
  float localEr;
  if (totalBases == 0) {
    localEr = 0;
  } else {
    localEr = ((float)numncMis)/totalBases;
  }
  
  cout << variant.chro << "\t" << variant.start << "\t" << variant.countAll << "\t" << variant.countPositive << "\t" << variant.countNegative << "\t" << variant.F1R2_all << "\t" << variant.F2R1_all << "\t" << variant.F1R2_alt << "\t" << variant.F2R1_alt << "\t" << variant.countAlt << "\t" << variant.countA << "\t" << variant.countAn << "\t" << variant.countC << "\t" << variant.countCn << "\t" << variant.countG << "\t" << variant.countGn << "\t" << variant.countT << "\t" << variant.countTn << "\t" << variant.inends << "\t" << variant.countJump << "\t" << setprecision(4) << fracBadMappingQual << "\t" << setprecision(2) << meanMis << "\t" << setprecision(2) << medianMis << "\t" << setprecision(3) << meanIndel << "\t" << setprecision(3) << medianIndel << "\t" << setprecision(3) << medianVRLength << "\t" << setprecision(2) << localEr << "\t" << variant.qualities << endl;

}

inline float CalcMedian (vector<unsigned int> &scores)
{
  float median;
  size_t size = scores.size();

  sort(scores.begin(), scores.end());

  if (size % 2 == 0)
  {
    median = ((float)scores[size / 2 - 1] + (float)scores[size / 2]) / 2;
  }
  else 
  {
    median = (float)scores[size / 2];
  }

  return median;
}
