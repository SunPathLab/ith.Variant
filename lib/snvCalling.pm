package snvCalling;

use strict;
use List::Util qw[min max sum];

#
# snv Calling
#


sub muTectCalling {

  my ($class, $muTectBin, $BAM, $NORMALBAM, $gfasta, $COSMIC, $DBSNP, $muTectOut, $vcfOut, $chrProcess, $mem, $tmpDir, $javaBin) = @_;
  my $cosmicOption = ($COSMIC eq '')? '':"--cosmic $COSMIC";

  my $cmd = "$javaBin -Xmx$mem -Djava.io.tmpdir=$tmpDir -jar $muTectBin -rf BadCigar --analysis_type MuTect --reference_sequence $gfasta $cosmicOption --dbsnp $DBSNP --validation_strictness LENIENT --input_file:normal $NORMALBAM --input_file:tumor $BAM --enable_extended_output --out $muTectOut -vcf $vcfOut";
  if ($chrProcess ne 'SRP') {
    $cmd = "$javaBin -Xmx$mem -Djava.io.tmpdir=$tmpDir -jar $muTectBin -rf BadCigar --analysis_type MuTect --reference_sequence $gfasta --intervals $chrProcess $cosmicOption --dbsnp $DBSNP --validation_strictness SILENT --input_file:normal $NORMALBAM --input_file:tumor $BAM --enable_extended_output --out $muTectOut -vcf $vcfOut";
  }

  return $cmd;

}


sub samtoolsCalling {

  my ($class, $samtoolsBin, $bcftoolsBin, $BAM, $NORMALBAM, $gfasta, $vcfOut, $maxDepth, $ignoreRG, $chrProcess, $sensitivity, $CfiftyOff) = @_;

  my $ignoreRGopt = '';
  if ($ignoreRG == 1) {
    $ignoreRGopt = '--ignore-RG';
  }

  my $regionOpt = '';
  if ($chrProcess ne 'SRP') {
    $regionOpt = '--region '.$chrProcess;
  }

  my $CfiftyOpt = '-C 50';
  if ($CfiftyOff == 1) {
    $CfiftyOpt = '';
  }

  my $sensOpt = '-P '.$sensitivity;

  my $cmd = "$samtoolsBin mpileup $ignoreRGopt $regionOpt -Eugd $maxDepth -t DP,SP -q 0 $CfiftyOpt -f $gfasta $BAM $NORMALBAM | $bcftoolsBin call -p 0.9 $sensOpt -vcf GQ - >$vcfOut";

  if ( $BAM eq $NORMALBAM ) {  #only one bam single calling
    $cmd = "$samtoolsBin mpileup $ignoreRGopt $regionOpt -Eugd $maxDepth -t DP,SP -q 0 $CfiftyOpt -f $gfasta $BAM | $bcftoolsBin call -p 0.9 $sensOpt -vcf GQ - >$vcfOut";
  }

  return $cmd;

}

sub strelkaCalling1 {

  my ($class, $strelkaBin, $NORMALBAM, $BAM, $GFASTA, $CONFIG, $outDir) = @_;

  my $cmd = "perl $strelkaBin --normal=$NORMALBAM --tumor=$BAM --ref=$GFASTA --config=$CONFIG --output-dir=$outDir";

  return $cmd;

}


sub strelkaCalling2 {

  my ($class, $outDir, $threads) = @_;

  my $cmd = "make -C $outDir -j $threads";

  return $cmd;

}


sub muTect2vcf {

  my ($class, $mutect2vcfBin, $inMutect, $outVCF, $samplePairNames) = @_;

  my $cmd = "perl $mutect2vcfBin $inMutect >$outVCF";

  if ($samplePairNames ne 'SRP') {       #redefined
    $cmd = "perl $mutect2vcfBin $inMutect $samplePairNames >$outVCF";
  }

  return $cmd;

}


sub vcfSort {

  my ($class, $vcfSortBin, $inVCF, $outVCF) = @_;

  my $cmd = "cat $inVCF | $vcfSortBin >$outVCF";

  return $cmd;

}


sub vcfFixGT {

  my ($class, $vcfFixGTbin, $inVCF) = @_;

  my $cmd = "perl $vcfFixGTbin $inVCF";

  return $cmd;

}


sub runAnnovar {

  my ($class, $annovarBin, $inVCF, $ANNOVARDB, $species) = @_;

  my $cmd = "perl $annovarBin $inVCF $ANNOVARDB -buildver $species -remove -protocol refGene,cytoBand,genomicSuperDups,phastConsElements46way,tfbsConsSites,gwasCatalog,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,exac03,popfreq_max_20150413,ljb26_all,clinvar_20150629 -operation g,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput --outfile $inVCF";

  if ($species eq 'hg38') {  #hg38, tfbsConsSites and popfreq_max not used
    $cmd = "perl $annovarBin $inVCF $ANNOVARDB -buildver $species -remove -protocol refGene,cytoBand,genomicSuperDups,phastConsElements100way,gwasCatalog,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,avsnp147,exac03,dbnsfp30a,clinvar_20150629 -operation g,r,r,r,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput --outfile $inVCF";
  }
  return $cmd;
}


sub convertVCFannovar {

  my ($class, $convertVCFannovarBin, $vcfMultiAnno, $vcfMultiAnnoVCF) = @_;

  my $cmd = "perl $convertVCFannovarBin $vcfMultiAnno $vcfMultiAnnoVCF";

  return $cmd;

}


sub grepSNVvcf {

  my ($class, $vcfMultiAnnoVCF, $vcfMultiAnnoVCFsnv, $pattern) = @_;

  my $cmd = "awk -F\"\\t\" \'\$1 \~ \/\^\#\/ \|\| \$8 \!\~ \/\^INDEL\/\' $vcfMultiAnnoVCF >$vcfMultiAnnoVCFsnv";
  if ($pattern ne 'SRP'){
    $cmd = "awk -F\"\\t\" \'\$1 \~ \/\^\#\/ \|\| \$8 \~ \/$pattern\/\' $vcfMultiAnnoVCF >$vcfMultiAnnoVCFsnv";
  }

  return $cmd;

}

sub grepINDELvcf {

  my ($class, $vcfMultiAnnoVCF, $vcfMultiAnnoVCFindel, $pattern) = @_;

  my $cmd = "awk -F\"\\t\" \'\$1 \~ \/\^\#\/ \|\| \$8 \~ \/\^INDEL\/\' $vcfMultiAnnoVCF >$vcfMultiAnnoVCFindel";
  if ($pattern ne 'SRP'){
    $cmd = "awk -F\"\\t\" \'\$1 \~ \/\^\#\/ \|\| \$8 \~ \/$pattern\/\' $vcfMultiAnnoVCF >$vcfMultiAnnoVCFindel";
  }

  return $cmd;

}

sub rechecksnv {

  my ($class, $rechecksnvBin, $recheckTable, $BAM, $recheckOut, $chrPref, $skipPileup) = @_;

  my $skipPileupOpt = ($skipPileup eq 'yes')? '--skipPileup' : '';
  my $cmd = "$rechecksnvBin --var $recheckTable --mapping $BAM $skipPileupOpt >$recheckOut";
  if ($chrPref ne 'SRP'){
    $cmd = "$rechecksnvBin --var $recheckTable --mapping $BAM $skipPileupOpt --chr $chrPref >$recheckOut";
  }

  return $cmd;

}

sub calTumorLOD {
  my ($class, $e, $le, $eg, $f, $v, $d) = @_;
  my $r = $d-$v;                         #reference allele count
  my @er = split("", $e);
  my @pms;
  my @pms0;
  foreach my $er (@er) {
    #print STDERR "$er\t";
    my $err = ord($er)-33;                     #ASCII -> quality
    $err -= 31 if ($err >= 64);                #qual shift
    $err = 62 if ($err > 62);                  #qual max
    $err = 10**(-$err/10);
    $er = max($err, $le);                      #maximum between local and phred score
    $er = 0.01 if ($er > 0.01);                #minimum if wierd case
    my $Pm = $f*(1-$er) + (1-$f)*($er/3);
    my $Pm0 = $er/3;
    #print STDERR "$er\t$Pm\t$Pm0\t";
    push(@pms, $Pm);
    push(@pms0, $Pm0);
  }
  my $Pr = $f*($eg/3) + (1-$f)*(1-$eg);
  my $lmut = ($Pr**$r)*(product(\@pms));

  my $Pr0 = 1-$eg;
  my $lref = ($Pr0**$r)*(product(\@pms0));

  if ($lref == 0) {
    return(100);
  } elsif ($lmut == 0) {
    return(-100);
  }
  my $lod = sprintf("%.6f", log10($lmut/$lref));
  if ($lod eq 'inf'){
    return(100);
  } elsif ($lod eq '-inf'){
    return(-100);
  }
  #print STDERR "\t$f\t$v\t$r\t$d\t$lmut\t$lref\t$lod\n";
  return($lod);
}

sub calNormalLOD {
  my ($class, $e, $le, $eg, $f, $v, $d) = @_;
  my $r = $d-$v;                        #reference allele count
  my $fg = 0.5;
  my @er = split("", $e);
  my @pms0;
  my @pmsg;
  foreach my $er (@er) {
    #print STDERR "$er\t";
    my $err = ord($er)-33;                     #ASCII -> quality
    $err -= 31 if ($err >= 64);                #qual shift
    $err = 62 if ($err > 62);                  #qual max
    $err = 10**(-$err/10);
    $er = max($err, $le);                      #maximum between local and phred score
    $er = 0.01 if ($er > 0.01);                #minimum if wierd case
    my $Pm0 = $er/3;
    my $Pmg = $fg*(1-$er) + (1-$fg)*($er/3);
    #print STDERR "$er\t$Pm0\t$Pmg\t";
    push(@pms0, $Pm0);
    push(@pmsg, $Pmg);
  }
  my $Pr0 = 1-$eg;
  my $lref = ($Pr0**$r)*(product(\@pms0));

  my $Prg = $fg*($eg/3) + (1-$fg)*(1-$eg);
  my $lger = ($Prg**$r)*(product(\@pmsg));

  if ($lref == 0) {
    return(-100);
  } elsif ($lger == 0) {
    return(100);
  }
  my $lod = sprintf("%.6f", log10($lref/$lger));
  if ($lod eq 'inf'){
    return(100);
  } elsif ($lod eq '-inf'){
    return(-100);
  }
  #print STDERR "\t$f\t$v\t$r\t$d\t$lref\t$lger\t$lod\n";
  return($lod);
}

sub log10 {
  my $n = shift;
  return log($n)/log(10);
}

sub product {
  my $a = shift;
  my $prod = 1;
  foreach my $element (@{$a}) {
    $prod *= $element;
  }
  return $prod;
}


1;


