use strict;
use Data::Dumper;
use File::Glob ':glob';
use File::Basename;

my $data = shift;
my $somaticInfo = shift;
my $genome = shift;
my $refcontextlength = 10;
my $split = 1;

my %somatic;
my %germline;  #may have multiple tumors
if ($somaticInfo ne '' and -s "$somaticInfo") {

  open IN, "$somaticInfo";
  while ( <IN> ){
    chomp;
    s/[\s\n]$//;
    my @columns = split /\t/;
    my $tumor = $columns[0];
    my $normal = $columns[1];

    $somatic{$tumor} = $normal;
    push(@{$germline{$normal}}, $tumor) if $normal ne 'undef';
  }
  close IN;
  print STDERR Dumper (\%somatic);
  print STDERR Dumper (\%germline);

}


open HS, "$genome";
my %genome = ();
my $chr = undef;
while ( <HS> ) {
  if (/^>(\w+).*?\n$/) {
    $chr = $1;
    $chr =~ s/^chr//;
  } else {
    s/\n//g;
    s/\s//g;
    $genome{$chr}.=$_;
  }
}
close HS;
print STDERR "genome loaded\n";


my $outdir = dirname($data)."/ToxoG/";
print STDERR "$outdir\n";
unless (-e "$outdir") {
  system("mkdir -p $outdir");
}


if ($split == 1) {

  my %colnames;
  my %colindex;
  my %fhs;
  open IN, "$data";
  while (<IN>) {
    chomp;
    my @cols = split /\t/;
    if ($_ =~ /^[\#]?chr\t/) {
      $_ =~ s/^\#//;
      for (my $i = 0; $i <= $#cols; $i++) {
        $colnames{$i} = $cols[$i];
        $colindex{$cols[$i]} = $i;
      }
      next;
    } else {
      my $chr = $cols[$colindex{'chr'}];
      my $pos = $cols[$colindex{'pos'}];
      my $ref = $cols[$colindex{'ref'}];
      my $alt = $cols[$colindex{'alt'}];
      my $refcontext = substr($genome{$chr}, ($pos-($refcontextlength+1)), 2*$refcontextlength+1);
      for (my $i = 0; $i <= $#cols; $i++) {
        if ($colnames{$i} =~ /^(.+?)maf$/) { #now it is sample maf

          my $sample = $1;
          my $normalsample = $somatic{$sample};

          if (exists($somatic{$sample})) { #it is a tumor sample

            if ($cols[$colindex{'somatic'}] =~ /$sample\[/) { #called somatic

              if ($cols[$i] =~ /\|/) { #split the var surrounding information
                my @infos = split(/\|/, $cols[$i]);
                my $maf = $infos[0];
                my $endsratio = $infos[1];
                my ($cmean, $cmedian) = split(',', $infos[2]);
                my @strandRatio = split(',', $infos[3]);
                my $strandRatio = $strandRatio[0];
                my $strandRatioRef = $strandRatio[1];
                my $strandFisherP = $strandRatio[2];
                my $badQualFrac = $infos[4];
                my $lod = $infos[5];
                my ($F1R2all, $F2R1all, $F1R2alt, $F2R1alt) = split(',', $infos[6]);

                #my $fh = $sample;
                unless (-e "$outdir/$sample\_toxog") {
                  open ( my $fh, ">>", "$outdir/$sample\_toxog" )  || die $!;
                  $fhs{$sample} = $fh;
                  print {$fhs{$sample}} "Chromosome\tStart_position\tEnd_position\tReference_Allele\tTumor_Seq_Allele2\ti_picard_oxoQ\tTumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tref_context\ti_t_ALT_F1R2\ti_t_ALT_F2R1\ti_t_REF_F1R2\ti_t_REF_F2R1\ti_t_Foxog\tVariant_Type\tTumor_Seq_Allele1\n";
                }
                my $F1R2ref = $F1R2all - $F1R2alt;
                my $F2R1ref = $F2R1all - $F2R1alt;
                my $Foxog = 0;
                if (($F2R1alt + $F1R2alt) > 0) {
                  $Foxog = ($ref =~ /[CA]/)? $F2R1alt/($F2R1alt + $F1R2alt) : $F1R2alt/($F2R1alt + $F1R2alt);
                }
                my $printout = join("\t", $chr, $pos, $pos, $ref, $alt, 0, $sample, $normalsample, $refcontext, $F1R2alt, $F2R1alt, $F1R2ref, $F2R1ref, $Foxog, 'SNP', $alt);
                print {$fhs{$sample}} "$printout\n";

              }                 #with info seperated by |
            }                   #called
          }                     #tumor
        }                       #maf
      }                         #each col
    }                           #each non header
  } #each line
  close IN;

} #split samples


sub round {
  my $number = shift;
  my $tmp = int($number);
  if ($number >= ($tmp+0.5)){
    $tmp++;
  }
  return $tmp;
}
