#!/usr/bin/env perl

#date   : 2025-11-21
#author : Xiangjian Gou

#load the modules
use strict;
use warnings;
use Getopt::Long;

#clear buffer
$| = 1;

#record version information
my $VERSION = 'CQD v1.0';

#set default options
my $inputGwasFile = "";
my $targetColIndex = "1,2,3,4"; #snp_id,chrom,pos,pvalue
my $mannerForColSep = "comma"; #comma, tab, space
my $outputFile = "candidate_qtls.txt";
my $pvalueThreshold = 0.05; #0.01, 0.05, 1
my $snpNumber = 1_000_000;
my $distanceBetweenSnp = 100_000;
my $clusterSnpNumber = 3;
my $annotationFile = "";
my $range = 100_000;
my $namePrefixOfChr = "";
my $nameSuffixOfChr = "";
my $version;
my $help;

#get options from command line
GetOptions(
    'inputGwasFile=s'       => \$inputGwasFile,
    'targetColIndex=s'      => \$targetColIndex,
    'mannerForColSep=s'     => \$mannerForColSep,
    'outputFile=s'          => \$outputFile,
    'pvalueThreshold=f'     => \$pvalueThreshold,
    'snpNumber=i'           => \$snpNumber,
    'distanceBetweenSnp=i'  => \$distanceBetweenSnp,
    'clusterSnpNumber=i'    => \$clusterSnpNumber,
    'annotationFile=s'      => \$annotationFile,
    'range=i'               => \$range,
    'namePrefixOfChr|np=s'  => \$namePrefixOfChr,
    'nameSuffixOfChr|ns=s'  => \$nameSuffixOfChr,
    'version+'              => \$version,
    'help+'                 => \$help,
);

#describe program information
my $usage = <<__GUIDE__;
##########################################################################################
Name:
  CQD - Candidate Qtl Detector (Detection of candidate QTLs based on the GWAS results)

Author:
  Xiangjian Gou (xjgou\@mail.hzau.edu.cn or 862137261\@qq.com)

Usage:
  perl CQD.pl option1 <value1> option2 <value2> ... optionN <valueN>

Options:
  -i  | -inputGwasFile   <STRING> : provide a GWAS file (must be provided !)
  -t  | -targetColIndex  <STRING> : column (snp_id,chrom,pos,pvalue) index in GWAS file (default: 1,2,3,4)
  -m  | -mannerForColSep <STRING> : delimiter (comma/tab/space) among column in GWAS file (default: comma)
  -o  | -outputFile      <STRING> : set the output file name (default: candidate_qtls.txt)
  -p  | -pvalueThreshold  <FLOAT> : p.value for Bonferroni correction (default: 0.05)
  -s  | -snpNumber          <INT> : number of SNPs for Bonferroni correction (default: 1000000)
  -d  | -distanceBetweenSnp <INT> : max distance between adjacent SNPs within a QTL (default: 100000)
  -c  | -clusterSnpNumber   <INT> : min number of significant SNPs within a QTL (default: 3)
  -a  | -annotationFile  <STRING> : a GFF3 file, if provided, will list candidate genes for each QTL
  -r  | -range              <INT> : up/downstream range of QTL to list candidate genes (default: 100000)
  -np | -namePrefixOfChr <STRING> : add a prefix to the chromosome in GWAS file (default: no prefix)
  -ns | -nameSuffixOfChr <STRING> : add a suffix to the chromosome in GWAS file (default: no suffix)
  -v  | -version                  : show the version information
  -h  | -help                     : show the help information

Note:
  1. Options '-t' and '-m' determine the format of the GWAS file.
  2. Options '-p' and '-s' determine the threshold (p/s) for identifying significant SNPs.
  3. Options '-d' and '-c' determine the definition criteria of a QTL.
  4. Options '-a' and '-r' enable the listing of candidate genes for a QTL.
  5. If '-r' is set to 0, candidate genes are only detected within the QTL.
  6. Options '-np' and '-ns' ensure that chromosome ID is consistent in the GWAS and GFF3 files.

Example:
  perl CQD.pl -i gwas_file.csv
##########################################################################################

__GUIDE__

#check the options
die "$VERSION\n" if $version;
die $usage if $help;
die "Error: option '-i' must be provided (plase check the help '-h') !\n" if ! $inputGwasFile;
die "Error: option '-t' has an improperly formatted assignment !\n" if $targetColIndex !~ /\d+,\d+,\d+,\d+/;
die "Error: option '-m' only be set to comma, tab, or space !\n" if $mannerForColSep ne 'comma' and $mannerForColSep ne 'tab' and $mannerForColSep ne 'space';
die "Error: option '-p' must be bigger than 0 !\n" if $pvalueThreshold <= 0;
die "Error: option '-s' must be an integer bigger than 0 !\n" if $snpNumber <= 0;
die "Error: option '-d' must be an integer bigger than 0 !\n" if $distanceBetweenSnp <= 0;
die "Error: option '-c' must be an integer bigger than 1 !\n" if $clusterSnpNumber <= 1;
die "Error: option '-r' cannot be set to a negative number !\n" if $range < 0;

#==============================start main program==============================

#set the significant threshold
my $bonferroniThreshold = $pvalueThreshold/$snpNumber;

#get the significant SNPs
my $significantSnps = getSigSnp($inputGwasFile, $bonferroniThreshold, $targetColIndex, $mannerForColSep);

#get the QTLs
my %qtls;
my $rank = 1;
foreach my $chr (sort keys %$significantSnps) {
    my @data = sort { $a->[0] <=> $b->[0] } @{$significantSnps->{$chr}};
    if (@data < $clusterSnpNumber) {
        next;
    }
    else {
        my %tmp;
        foreach my $i (0 .. $#data-1) {
            my ($pos1, $snpID1, $pvalue1) = @{$data[$i]};
            my ($pos2, $snpID2, $pvalue2) = @{$data[$i+1]};
            my $distance = abs($pos2-$pos1);
            if ($distance <= $distanceBetweenSnp) {
                $tmp{$snpID1} = [$pos1, $pvalue1];
                $tmp{$snpID2} = [$pos2, $pvalue2];
                #deal with the last one
                if ($i == $#data-1) {
                    my $snpNum = keys %tmp;
                    if ($snpNum >= $clusterSnpNumber) {
                        my @qtlInfo = getQtlInfo(\%tmp);
                        $qtls{$rank} = [$chr, $snpNum, @qtlInfo];
                        $rank++;
                    }
                    undef %tmp;
                }
            }
            else {
                my $snpNum = keys %tmp;
                if ($snpNum >= $clusterSnpNumber) {
                    my @qtlInfo = getQtlInfo(\%tmp);
                    $qtls{$rank} = [$chr, $snpNum, @qtlInfo];
                    $rank++;
                }
                undef %tmp;
            }
        }
    }
}

#output the QTLs
if ($annotationFile eq '') {
    open my $out, '>', $outputFile;
    print $out join("\t", qw/chr lead_snp_id lead_snp_pos lead_snp_pvalue qtl_start qtl_end snp_number_in_qtl/), "\n";
    foreach my $rank (sort { $qtls{$a}[4] <=> $qtls{$b}[4] } keys %qtls) { #sorted by pvalue
        my ($chr, $snpNum, $leadSnpID, $leadSnpPos, $leadSnpPvalue, $minPos, $maxPos) = @{$qtls{$rank}};
        print $out join("\t", $chr, $leadSnpID, $leadSnpPos, $leadSnpPvalue, $minPos, $maxPos, $snpNum), "\n";
    }
    close $out;
}
else {
    my $geneInfo = getGeneInfo($annotationFile);
    open my $out, '>', $outputFile;
    print $out join("\t", qw/chr lead_snp_id lead_snp_pos lead_snp_pvalue qtl_start qtl_end snp_number_in_qtl candidate_gene_number candidate_genes/), "\n";
    foreach my $rank (sort { $qtls{$a}[4] <=> $qtls{$b}[4] } keys %qtls) { #sorted by pvalue
        my ($chr, $snpNum, $leadSnpID, $leadSnpPos, $leadSnpPvalue, $minPos, $maxPos) = @{$qtls{$rank}};
        my $chrInGff3 = $namePrefixOfChr.$chr.$nameSuffixOfChr;
        my $candidates = getAImGenes($minPos, $maxPos, $range, $geneInfo->{$chrInGff3});
        my $geneNumber = @$candidates;
        my $geneSet = $geneNumber == 0 ? 'NA' : join(";", @$candidates);
        print $out join("\t", $chr, $leadSnpID, $leadSnpPos, $leadSnpPvalue, $minPos, $maxPos, $snpNum, $geneNumber, $geneSet), "\n";
    }
    close $out;
}

#================================================================================

sub getSigSnp {
    my ($file, $bonferroniThreshold, $targetColIndex, $mannerForColSep) = @_;

    #parse the file format
    my ($index1, $index2, $index3, $index4) = split /,/, $targetColIndex;
    $index1 -= 1; #snp_id
    $index2 -= 1; #chrom
    $index3 -= 1; #pos
    $index4 -= 1; #pvalue
    my $splitChar = do {
        if    ($mannerForColSep eq 'comma') { ","   }
        elsif ($mannerForColSep eq 'tab'  ) { "\\t" }
        elsif ($mannerForColSep eq 'space') { " "   }
        else  { die "Error: bug1 at code !\n"       }
    };

    #determine if the GWAS file has been sorted according to the pvalue
    open my $in1, '<', $file or die "Error: cannot open file '$file': $!";
    <$in1>;
    my $count = 0;
    my $firstLineNumber = 1000;
    my @pvalues;
    while (<$in1>) {
        $count++;
        s/[\r\n]+//;
        my @row = split /$splitChar/;
        die "Error: snp_id index out of range, please check the options '-t' and '-m' !\n" if $index1 > $#row;
        die "Error: chrom index out of range, please check the options '-t' and '-m' !\n"  if $index2 > $#row;
        die "Error: pos index out of range, please check the options '-t' and '-m' !\n"    if $index3 > $#row;
        die "Error: pvalue index out of range, please check the options '-t' and '-m' !\n" if $index4 > $#row;
        my $pvalue = $row[$index4];
        $pvalue = 1 if $pvalue =~ /NA/i;
        die "Error: pvalue index is not right !\n" if $pvalue > 1 or $pvalue < 0;
        $count <= $firstLineNumber ? push(@pvalues, $pvalue) : last;
    }
    close $in1;
    my $judge = "yes";
    if (@pvalues == 0) {
        die "Error: the GWAS file has no SNP !\n";
    }
    elsif (@pvalues == 1) {
        die "Error: the GWAS file has 1 SNP !\n";
    }
    else {
        foreach my $i (0 .. $#pvalues-1) {
            if ($pvalues[$i] > $pvalues[$i+1]) {
                $judge = "no";
                last;
            }
        }
    }
    print "LOG: GWAS file has been sorted according to the pvalue: $judge\n";

    #get the significant SNPs
    open my $in2, '<', $file or die "Error: cannot open file '$file': $!";
    <$in2>;
    my %significantSnps;
    while (<$in2>) {
        s/[\r\n]+//;
        my @row = split /$splitChar/;
        my ($snpID, $chr, $pos, $pvalue) = @row[$index1, $index2, $index3, $index4];
        $pvalue = 1 if $pvalue =~ /NA/i;
        if ($pvalue <= $bonferroniThreshold) {
            push @{$significantSnps{$chr}}, [$pos, $snpID, $pvalue];
        }
        else {
            last if $judge eq 'yes';
        }
    }
    close $in2;
    print "Warning: there is no significant SNP !\n" if keys %significantSnps == 0;

    #output
    return \%significantSnps;
}

sub getQtlInfo {
    my $data = shift;
    my @snpIDs = sort keys %$data;
    my $leadSnpID     = $snpIDs[0];
    my $leadSnpPos    = $data->{$snpIDs[0]}[0];
    my $leadSnpPvalue = $data->{$snpIDs[0]}[1];
    my $minPos = $data->{$snpIDs[0]}[0];
    my $maxPos = $data->{$snpIDs[0]}[0];
    foreach my $i (1 .. $#snpIDs) {
        my $snpID  = $snpIDs[$i];
        my $pos    = $data->{$snpID}[0];
        my $pvalue = $data->{$snpID}[1];
        $leadSnpID     = $snpID  if $leadSnpPvalue > $pvalue;
        $leadSnpPos    = $pos    if $leadSnpPvalue > $pvalue;
        $leadSnpPvalue = $pvalue if $leadSnpPvalue > $pvalue;
        $minPos = $pos if $minPos > $pos;
        $maxPos = $pos if $maxPos < $pos;
    }
    return $leadSnpID, $leadSnpPos, $leadSnpPvalue, $minPos, $maxPos;
}

sub getGeneInfo {
    my $file = shift;
    open my $in, '<', $file or die "Error: cannot open file '$file': $!";
    my %genes;
    while (<$in>) {
        s/[\r\n]+//;
        next if /\A#/;
        my ($chr, $name, $start, $end, $info) = (split /\t/, $_)[0, 2, 3, 4, 8];
        next if $name !~ /gene/;
        my ($geneID) = $info =~ /\AID=gene:([^;]+);/;
        $genes{$chr}{$geneID} = [$start, $end];
    }
    close $in;
    return \%genes;
}

sub getAImGenes {
    my ($qtlStart, $qtlEnd, $range, $genes) = @_;
    my $qtlStart2 = $qtlStart-$range;
    my $qtlEnd2 = $qtlEnd+$range;
    my @output;
    foreach my $geneID (sort { $genes->{$a}[0] <=> $genes->{$b}[0] } keys %$genes) {
        my ($geneStart, $geneEnd) = @{$genes->{$geneID}};
        if ($geneStart > $qtlEnd2) {
            last;
        }
        elsif ($geneEnd < $qtlStart2) {
            next;
        }
        else {
            push @output, $geneID;
        }
    }
    return \@output;
}
