#!/usr/bin/perl
use strict;
use File::Basename;

my $fin = shift;
my $vcfOriginal = shift;
(my $vcfOut = $vcfOriginal) =~ s/\.vcf$/\.mod\.vcf/;

if ($fin eq '' or $vcfOriginal eq '') {
  print STDERR "the annovar table and vcf must be specified!!\n";
  exit 22;
}

open(ANNOVARTABLE, "<" . $fin) || die("Could not open $fin file!");
open(VCFORI, "<" . $vcfOriginal) || die("Could not open original VCF file!");
open(OUT, ">" .  $vcfOut)|| die("Could not open $vcfOut  file!");

my $headerAdd = '';
$headerAdd .= '##INFO=<ID=ESP6500,Number=1,Type=Integer,Description="ESP6500 score as reported by Annovar">'."\n";
$headerAdd .= '##INFO=<ID=1KG,Number=1,Type=Integer,Description="1KG Oct2014 score as reported by Annovar">'."\n";
$headerAdd .= '##INFO=<ID=1KG.AFR,Number=1,Type=Integer,Description="1KG African Oct2014 score as reported by Annovar">'."\n";
$headerAdd .= '##INFO=<ID=1KG.EAS,Number=1,Type=Integer,Description="1KG East Asian Oct2014 score as reported by Annovar">'."\n";
$headerAdd .= '##INFO=<ID=1KG.EUR,Number=1,Type=Integer,Description="1KG European Oct2014 score as reported by Annovar">'."\n";
$headerAdd .= '##INFO=<ID=dbSNP,Number=0,Type=Flag,Description="dbSNP138 membership as reported by Annovar">'."\n";
$headerAdd .= '##INFO=<ID=ExAC_ALL,Number=.,Type=String,Description="ExAC_ALL annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=ExAC_AFR,Number=.,Type=String,Description="ExAC_AFR annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=ExAC_AMR,Number=.,Type=String,Description="ExAC_AMR annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=ExAC_EAS,Number=.,Type=String,Description="ExAC_EAS annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=ExAC_FIN,Number=.,Type=String,Description="ExAC_FIN annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=ExAC_NFE,Number=.,Type=String,Description="ExAC_NFE annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=ExAC_OTH,Number=.,Type=String,Description="ExAC_OTH annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=ExAC_SAS,Number=.,Type=String,Description="ExAC_SAS annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=PopFreqMax,Number=.,Type=String,Description="PopFreqMax annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=function,Number=1,Type=String,Description="function as reported by Annovar">'."\n";
$headerAdd .= '##INFO=<ID=functionClass,Number=1,Type=String,Description="functional Class as reported by Annovar">'."\n";
$headerAdd .= '##INFO=<ID=geneName,Number=1,Type=String,Description="gene name as reported by Annovar">'."\n";
$headerAdd .= '##INFO=<ID=AAChange,Number=1,Type=String,Description="Amino Acid change as reported by Annovar">'."\n";
$headerAdd .= '##INFO=<ID=segdup.score,Number=1,Type=Integer,Description="Segdup score as reported by Annovar">'."\n";
$headerAdd .= '##INFO=<ID=cytoBand,Number=.,Type=String,Description="cytoBand annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=phastCons.score,Number=.,Type=String,Description="phastConsElements46way annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=phastCons.lod,Number=.,Type=String,Description="phastConsElements46way annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=tfbsConsSites.score,Number=1,Type=Integer,Description="tfbsConsSites annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=tfbsConsSites.name,Number=.,Type=String,Description="tfbsConsSites annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=dgvMerged,Number=.,Type=String,Description="dgvMerged annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=gwasCatalog,Number=.,Type=String,Description="gwasCatalog annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=SIFT_score,Number=.,Type=String,Description="SIFT_score annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=SIFT_pred,Number=.,Type=String,Description="SIFT_pred annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=Polyphen2_HDIV_score,Number=.,Type=String,Description="Polyphen2_HDIV_score annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=Polyphen2_HDIV_pred,Number=.,Type=String,Description="Polyphen2_HDIV_pred annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=Polyphen2_HVAR_score,Number=.,Type=String,Description="Polyphen2_HVAR_score annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=Polyphen2_HVAR_pred,Number=.,Type=String,Description="Polyphen2_HVAR_pred annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=LRT_score,Number=.,Type=String,Description="LRT_score annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=LRT_pred,Number=.,Type=String,Description="LRT_pred annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=MutationTaster_score,Number=.,Type=String,Description="MutationTaster_score annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=MutationTaster_pred,Number=.,Type=String,Description="MutationTaster_pred annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=MutationAssessor_score,Number=.,Type=String,Description="MutationAssessor_score annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=MutationAssessor_pred,Number=.,Type=String,Description="MutationAssessor_pred annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=FATHMM_score,Number=.,Type=String,Description="FATHMM_score annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=FATHMM_pred,Number=.,Type=String,Description="FATHMM_pred annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=RadialSVM_score,Number=.,Type=String,Description="RadialSVM_score annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=RadialSVM_pred,Number=.,Type=String,Description="RadialSVM_pred annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=LR_score,Number=.,Type=String,Description="LR_score annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=LR_pred,Number=.,Type=String,Description="LR_pred annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=VEST3_score,Number=.,Type=String,Description="VEST3_score annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=CADD_raw,Number=.,Type=String,Description="CADD_raw annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=CADD_phred,Number=.,Type=String,Description="CADD_phred annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=GERP++_RS,Number=.,Type=String,Description="GERP++_RS annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=phyloP46way_placental,Number=.,Type=String,Description="phyloP46way_placental annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=phyloP100way_vertebrate,Number=.,Type=String,Description="phyloP100way_vertebrate annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=SiPhy_29way_logOdds,Number=.,Type=String,Description="SiPhy_29way_logOdds annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=clinvar_20150629,Number=.,Type=String,Description="clinvar_20150629 annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=clinvar_20150330,Number=.,Type=String,Description="clinvar_20150330 annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=ALLELE_END,Number=0,Type=Flag,Description="Flag the end of ANNOVAR annotation for one alternative allele">'."\n";
$headerAdd .= '##INFO=<ID=nci60,Number=1,Type=Float,Description="nci60 annotation provided by ANNOVAR">'."\n";
$headerAdd .= '##INFO=<ID=avsnp144,Number=.,Type=String,Description="avsnp144 annotation provided by ANNOVAR">'."\n";


my %mapping;
$mapping{'Func.refGene'} = 'function';
$mapping{'Gene.refGene'} = 'geneName';
$mapping{'GeneDetail.refGene'} = 'geneDetail';
$mapping{'ExonicFunc.refGene'} = 'functionalClass';
$mapping{'AAChange.refGene'} = 'AAChange';
$mapping{'genomicSuperDups'} = 'segdup.score';
$mapping{'phastConsElements46way'} = 'phastCons'; #score and lod
$mapping{'phastConsElements100way'} = 'phastCons'; #score and lod
$mapping{'esp6500siv2_all'} = 'ESP6500';
$mapping{'1000g2015aug_all'} = '1KG';
$mapping{'1000g2015aug_all'} = '1KG';
$mapping{'1000g2014oct_afr'} = '1KG.AFR';
$mapping{'1000g2015aug_afr'} = '1KG.AFR';
$mapping{'1000g2014oct_eas'} = '1KG.EAS';
$mapping{'1000g2015aug_eas'} = '1KG.EAS';
$mapping{'1000g2014oct_eur'} = '1KG.EUR';
$mapping{'1000g2015aug_eur'} = '1KG.EUR';
$mapping{'1000g2014oct_eur'} = '1KG.SAS';
$mapping{'1000g2015aug_eur'} = '1KG.SAS';
$mapping{'snp138'} = 'dbSNP';
$mapping{'avsnp147'} = 'dbSNP';


my $headerbuffer;
my $headerColumns;
my $headerbufferFlag = 0;
while ( <VCFORI> ) {
  chomp;
  if ($_ =~ /^\#/){  #header
    if ( $_ =~ /^[#]+CHROM\t/ ) {
      $headerColumns = "$_\n";
      last;
    } else {
      if ( $_ =~ /^\#\#INFO\=\<ID\=ANNOVAR\_DATE/ ) {
        $headerbufferFlag = 1;
      }
      $headerbuffer .= "$_\n" if ($headerbufferFlag == 0);
    }
  }
}
close VCFORI;

print OUT "$headerbuffer";              #header part 1
print OUT "$headerAdd";                 #header part 2
print OUT "$headerColumns";             #header part 3


my %colindex;
my %colnames;
while ( <ANNOVARTABLE> ) {
  chomp;
  if ($_ =~ /^Chr\tStart/) {   #now it is the header
    my @cols = split /\t/;
    for(my $i = 0; $i <= $#cols; $i++) {
      $colindex{$cols[$i]} = $i;
      $colnames{$i} = $cols[$i];
    }
  } else {

    my @item = split /\t/;

    my $buffer;
    my @subItem;  #original vcf info
    for (my $j = $colindex{'Otherinfo'}+3; $j <= $#item; $j++) {
      push(@subItem, $item[$j]);
    }
    $buffer.= "$subItem[0]\t$subItem[1]\t";                                                                                                                        #Chr,Pos
    $buffer.= (exists($colindex{'snp138'}))? "$item[$colindex{'snp138'}]\t" : ((exists($colindex{'avsnp147'}))? "$item[$colindex{'avsnp147'}]\t": "\.\t");         #ID
    $buffer.="$item[$colindex{'Ref'}]\t$item[$colindex{'Alt'}]\t$subItem[5]\t$subItem[6]\t";                                                                       #Ref,Alt,Qual,Filter

    ## Update INFO with annotations
    my $info="$subItem[7]";
    $info =~ s/(\;dbSNP)?\;function\=.+?$//;   #erase the previous annotation
    $info .= ';';

    for (my $i = 0; $i <= $colindex{'Otherinfo'}; $i++) {  #each annotation by annovar
      next if ($colnames{$i} eq 'Chr');
      next if ($colnames{$i} eq 'Start');
      next if ($colnames{$i} eq 'End');
      next if ($colnames{$i} eq 'Ref');
      next if ($colnames{$i} eq 'Alt');
      next if ($colnames{$i} eq 'Otherinfo');
      next if ($colnames{$i} eq 'dgvMerged');

      next if ($item[$i] eq '.');    #skip no info line
      my $attr = (exists($mapping{$colnames{$i}}))? $mapping{$colnames{$i}}:$colnames{$i};

      #cols needs care
      if ($attr eq 'geneDetail') {
        $item[$i] =~ s/\;/\,/;
        if ($item[$i] =~ /^(dist\=.+?)\,(dist\=.+?)$/) {
          my $dist1 = $1;
          my $dist2 = $2;
          if ($info =~ /\,/) {
            $info =~ s/geneName\=(.+?)\,(.+?)\;/geneName\=\1\($dist1\)\,\2\($dist2\)\;/;
          } else {
            print STDERR "inconsistent geneName with dists: $item[$i]\t$info\n";
          }
          next;
        }
        $info .= "$attr\=$item[$i]\;";
      } elsif ($attr eq 'segdup.score'){
        $item[$i] =~ /Score\=(.+?)\;/;
        my $segdupscore = $1;
        $info .= "$attr\=$segdupscore\;";
      } elsif ($attr eq 'phastCons'){
        $item[$i] =~ /Score=(\d+)\;Name\=lod\=(\d+)/;
        my $pcscore = $1;
        my $pclod = $2;
        $info .= "$attr\.score\=$pcscore\;";
        $info .= "$attr\.lod\=$pclod\;";
      } elsif ($attr eq 'tfbsConsSites'){
        $item[$i] =~ /^Score=(\d+)\;Name\=(.+?)$/;
        my $tfscore = $1;
        my $tfname = $2;
        $info .= "$attr\.score\=$tfscore\;";
        $info .= "$attr\.name\=$tfname\;";
      } elsif ($attr eq 'gwasCatalog'){
        $item[$i] =~ /^Name=(.+?)$/;
        my $gwasname = $1;
        $info .= "$attr\=$gwasname\;";
      } elsif ($attr eq 'dbSNP'){
        $info .= "$attr\;";
      } elsif ($attr eq 'clinvar_20150629'){
        $info .= "$attr\;";
      }
      #cols NO care
      else {
        $info .= "$attr\=$item[$i]\;";
      }
    } #each annotation

    $buffer .= $info;

    ## Format & Data
    for (my $j = 8; $j <= $#subItem; $j++) {
      $buffer .= "\t".$subItem[$j];
    }
    print OUT "$buffer\n";
  }
}
close ANNOVARTABLE;
close OUT;


exit 0;
