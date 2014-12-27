#
use strict;

my ($targetFile,$flank) = @ARGV;
if(!defined $flank){
    $flank = 500;
}

open(TF, "<$targetFile") or die "Failed to read $targetFile\n";
open(GRP, ">$targetFile.GerpMatrix") or die "Failed to write $targetFile.GerpMatrix\n";
open(TRI, ">$targetFile.TriMatrix") or die "Failed to write $targetFile.GerpMatrix\n";

my $tabix = "tabix -s 1 -b 3 -e 3";
my @vecTemplate = split //, "0"x(($flank*2)+1);
while(<TF>){
    chomp;
    my @sp = split /\s+/;
    $sp[0] =~ s/chr//;
    my $leftC = $sp[1] - $flank;
    my $rightC = $sp[1] + $flank;
    my $slice = "$sp[0]:$leftC-$rightC";
    my @rv = `$tabix chr$sp[0].gerpVec.bed.filteredMap100.whiteList.bed.bgz $slice`;
    my @vec = @vecTemplate;
    my @tri = @vecTemplate;
    my $absScore = 0;
    for my $i (0 .. $#rv){
	my @rvsp = split /\s+/, $rv[$i];
	my $rp = $rvsp[2]-$leftC;
	$vec[$rp] = $rvsp[5];
	$tri[$rp] = $rvsp[4];
	$absScore += abs($rvsp[5]);
	#print "$rv[$i]";
    }
    next if ($absScore == 0);
    my $ov = join "\t", @vec;
    my $ot = join "\t", @tri;
    print GRP "$sp[0]\t$sp[1]\t$sp[3]\t$sp[4]\t$ov\n";
    print TRI "$sp[0]\t$sp[1]\t$sp[3]\t$sp[4]\t$ot\n";
}
