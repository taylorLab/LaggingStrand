# For every nucleotide, classify it into one of 64 groups based on
# flanking nucleotides.
# Martin Taylor 2014

use strict;

my ($fa) = @ARGV;

my ($tripos) = classes();
open (FA,"<$fa") or die "Failed to read $fa\n";




my ($name,$nom);
my $sq;
my $process=0;
while (<FA>){
    chomp;
    if(/^>(\S+)/){
	$nom = $1;
	$process=1;
    }
    else {
	$sq .= $_;
    }
    if(eof){
	$process=1;
    }
    if($process){
	if(defined($sq)){
	    $sq = uc $sq;
	    my ($classed) = triNuc($tripos,$sq);
	    for my $y (0 .. $#{$classed}){
		my $gid = ($y+1);
		print "$nom\t$y\t$classed->[$y]\n";
	    }
	}
	$name=$nom;
	$process=0;
	$sq = "";
    }
}

sub classes {
    my %cls;
    my @nucs = qw{A C G T};
    my $id=1;
    foreach my $n (@nucs){
	foreach my $m (@nucs){
	    foreach my $o (@nucs){
		my $tri = $n.$m.$o;
		if (!defined $cls{$tri}){
		    $cls{$tri} = $id++;
		}
	    }
	}
    }
    return(\%cls);
}

sub triNuc {
    my ($tref,$sqstream) = @_;
    my $ln = length($sqstream);
    my @catVec = split '', "0"x$ln;
    
    for my $i (0 .. ($ln-3)){
	my $trip=substr($sqstream,$i,3);
	if (!$tref->{$trip}){
	    #print STDERR "WARN: $trip is not defined\n";
	    
	}
	$catVec[$i+1]=$tref->{$trip};
    }
    return(\@catVec);
}

