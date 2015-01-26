#
# Produce a GID indexed genome array for the 5' end
# coordinate of bed annotation, counting ends and
# filtering on score.
#
# Martin Taylor, July 2013.

use strict;

my $scoreThresh = 20;
my @sacCer3ChromLengths = ( 0,
			    230218,
			    813184,
			    316620,
			    1531933,
			    576874,
			    270161,
			    1090940,
			    562643,
			    439888,
			    745751,
			    666816,
			    1078177,
			    924431,
			    784333,
			    1091291,
			    948066,
			    85779
    );


my %chrTranslate = ( 'chrI' => 1,
		     'chr1' => 1,
		     '1' => 1,
                     'chrII' => 2,
		     'chr2' => 2,
		     '2' => 2,
                     'chrIII' => 3,
		     'chr3' => 3,
		     '3' => 3,
                     'chrIV' => 4,
		     'chr4' => 4,
		     '4' => 4,
                     'chrV' => 5,
		     'chr5' => 5,
		     '5' => 5,
                     'chrVI' => 6,
		     'chr6' => 6,
		     '6' => 6,
                     'chrVII' => 7,
		     'chr7' => 7,
		     '7' => 7,
                     'chrVIII' => 8,
		     'chr8' => 8,
		     '8' => 8,
                     'chrIX' => 9,
		     'chr9' => 9,
		     '9' => 9,
                     'chrX' => 10,
		     'chr10' => 10,
		     '10' => 10,
                     'chrXI' => 11,
		     'chr11' => 11,
		     '11' => 11,
                     'chrXII' => 12,
		     'chr12' => 12,
		     '12' => 12,
                     'chrXIII' => 13,
		     'chr13' => 13,
		     '13' => 13,
                     'chrXIV' => 14,
		     'chr14' => 14,
		     '14' => 14,
                     'chrXV' => 15,
		     'chr15' => 15,
		     '15' => 15,
                     'chrXVI' => 16,
		     'chr16' => 16,
		     '16' => 16,
                     'chrM' => 17,
		     'chr17' => 17,
		     '17' => 17
    );

my @sacCer3ChromOffsets = (0,0);
for my $i (1 .. $#sacCer3ChromLengths){
    $sacCer3ChromOffsets[$i+1] = $sacCer3ChromOffsets[$i] + $sacCer3ChromLengths[$i];
}

my @gidPosFor = ();
my @gidPosRev = ();
my $vectorLength = $sacCer3ChromOffsets[$#sacCer3ChromOffsets];


for my $i (1 .. $#sacCer3ChromLengths){
    for my $j (1 .. $sacCer3ChromLengths[$i]){
	my $gNucPos = $sacCer3ChromOffsets[$i] + $j;
	$gidPosFor[$gNucPos] = 0;
    }
}

@gidPosRev = @gidPosFor;

while(<>){
    chomp;
    my @sp = split /\s+/;
    next if ($#sp < 2);
    next if ($sp[0] =~ /#/);
    next if ($sp[4] < $scoreThresh);
    if (!$chrTranslate{$sp[0]}){
	print STDERR "WARNING: $sp[0] not a recognised chromosome\n";
	next;
    }
    ###
    # With the riboseq protocol, the incorporated ribo is one nucelotide
    # upstream (5') of the 5' end of the read and on the opposite strand.
    # Output coordianates are GID (Genome InDex) 1-based.
    # As BED file input is zero-based, half open coordinates, we don't
    # need to subtract 1 from the +ve strand match, it's already correct
    # in 1-based coordinates.
    ###
    if ($sp[5] eq "+"){
	$gidPosRev[$sacCer3ChromOffsets[$chrTranslate{$sp[0]}]+($sp[1])]++;
    }
    else {
	$gidPosFor[$sacCer3ChromOffsets[$chrTranslate{$sp[0]}]+($sp[2]+1)]++;
    }
}

for my $i (1 .. $vectorLength){
    print "$i\t$gidPosFor[$i]\t$gidPosRev[$i]\n";
}

