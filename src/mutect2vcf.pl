use strict;

my $mutect = shift;
my $realSamples = shift;

my $vcfheader .= '##fileformat=VCFv4.1'."\n";
$vcfheader .= '##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">'."\n";
$vcfheader .= '##FILTER=<ID=REJECT,Description="Rejected as a confident somatic mutation">'."\n";
$vcfheader .= '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'."\n";
$vcfheader .= '##FORMAT=<ID=BQ,Number=A,Type=Float,Description="Average base quality for reads supporting alleles">'."\n";
$vcfheader .= '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">'."\n";
$vcfheader .= '##FORMAT=<ID=FA,Number=A,Type=Float,Description="Allele fraction of the alternate allele with regard to reference">'."\n";
$vcfheader .= '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'."\n";
$vcfheader .= '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'."\n";
$vcfheader .= '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">'."\n";
$vcfheader .= '##FORMAT=<ID=SS,Number=1,Type=Integer,Description="Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown">'."\n";
$vcfheader .= '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">'."\n";
$vcfheader .= '##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">'."\n";
$vcfheader .= '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic event">'."\n";
$vcfheader .= '##INFO=<ID=VT,Number=1,Type=String,Description="Variant type, can be SNP, INS or DEL">'."\n";

my $samples = 'SRP';
my %colindex;
if ($mutect =~ /\.gz$/) {
  open MU, "gzip -dc $mutect |";
} else {
  open MU, "$mutect";
}
while ( <MU> ) {
  chomp;
  next if /^#/;   #skip comments
  my @cols = split /\t/;
  if ( /^contig\t/ ){
    for (my $i = 0; $i <= $#cols; $i++) {
      $colindex{$cols[$i]} = $i;
    }
    next;
  }
  #my ($contig, $position, $context, $ref_allele, $alt_allele, $tumor_name, $normal_name, $score, $dbsnp_site, $covered, $power, $tumor_power, $normal_power, $normal_power_nsp, $normal_power_wsp, $total_pairs, $improper_pairs, $map_Q0_reads, $init_t_lod, $t_lod_fstar, $t_lod_lqs, $t_lod_fstar_forward, $t_lod_fstar_reverse, $tumor_f, $tumor_f_lb, $contaminant_fraction, $contaminant_lod, $t_q20_count, $t_ref_count, $t_alt_count, $t_ref_sum, $t_alt_sum, $t_ref_max_mapq, $t_alt_max_mapq, $t_ins_count, $t_del_count, $normal_best_gt, $init_n_lod, $n_lod_fstar, $normal_f, $normal_f_quals, $normal_artifact_lod_tf, $normal_artifact_lod_low_tf, $normal_artifact_lod_nf, $normal_artifact_lod_nfq, $n_q20_count, $n_ref_count, $n_alt_count, $n_ref_sum, $n_alt_sum, $power_to_detect_positive_strand_artifact, $power_to_detect_negative_strand_artifact, $strand_bias_counts, $tumor_alt_fpir_median, $tumor_alt_fpir_mad, $tumor_alt_rpir_median, $tumor_alt_rpir_mad, $alt_fpir, $alt_rpir, $powered_filters, $normal_artifact_power_tf, $normal_artifact_power_low_tf, $normal_artifact_power_nf, $normal_global_qll, $normal_local_qll, $normal_qmodel_lod, $observed_in_normals_count, $failure_reasons, $judgement) = split /\t/;

  if ($samples eq 'SRP'){
    $vcfheader .= "\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    my $tumor_name = $cols[$colindex{'tumor_name'}];
    my $normal_name = $cols[$colindex{'normal_name'}];
    if ($realSamples ne '') {
      my @realSamples = split (/\,/, $realSamples);
      $samples = join("\t", @realSamples);
      $vcfheader .= "\t$samples\n";
    } elsif ($tumor_name ne '' and $normal_name ne ''){
      $vcfheader .= "\t$tumor_name\t$normal_name\n";
      $samples = "$tumor_name\t$normal_name";
    }
    print "$vcfheader";
  }

  my $chr = $cols[$colindex{'contig'}];
  my $pos = $cols[$colindex{'position'}];
  my $id = '.';
  my $ref = $cols[$colindex{'ref_allele'}];
  my $alt = $cols[$colindex{'alt_allele'}];
  my $qual = '.';
  my $judgement = $cols[$colindex{'judgement'}];
  my $filter = ($judgement eq 'KEEP')? 'PASS':$judgement;
  my $t_lod_fstar = $cols[$colindex{'t_lod_fstar'}];
  my $init_n_lod = $cols[$colindex{'init_n_lod'}];
  my $info .= 't_lod_fstar='.$t_lod_fstar.';';
  $info .= 'init_n_lod='.$init_n_lod.';';
  my $format = 'GT:AD:BQ:DP:FA';

  my $t_ref_count = $cols[$colindex{'t_ref_count'}];
  my $t_alt_count = $cols[$colindex{'t_alt_count'}];
  my $n_ref_count = $cols[$colindex{'n_ref_count'}];
  my $n_alt_count = $cols[$colindex{'n_alt_count'}];
  my $t_alt_max_mapq = $cols[$colindex{'t_alt_max_mapq'}];
  my $t_dp = $t_ref_count + $t_alt_count;
  my $n_dp = $n_ref_count + $n_alt_count;
  my $tumor_f = $cols[$colindex{'tumor_f'}];
  my $normal_f = (exists($colindex{'normal_f'}))? $cols[$colindex{'normal_f'}]: (($n_dp > 0)? sprintf("%.5f", $n_alt_count/$n_dp):0);
  my $tumor = '0/1:'."$t_ref_count,$t_alt_count\:"."$t_alt_max_mapq\:"."$t_dp\:"."$tumor_f";
  my $normal = '0:'."$n_ref_count,$n_alt_count\:".'.:'."$n_dp\:"."$normal_f";
  printf("%s\n", join("\t", $chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $tumor, $normal));

}
close MU;

exit 0;
