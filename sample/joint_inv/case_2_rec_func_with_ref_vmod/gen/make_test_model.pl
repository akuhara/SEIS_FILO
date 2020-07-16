#!/usr/bin/perl
use strict;
use warnings;


my $ocean_thick = 2.0;
my $dz = 0.1;
my $z_max = 23.0;

# added anomaly 
my @zs  = (6, 12, 18);
my @dvs = (-0.05, -0.2, 0.1, 0.05); 

# Refernece model 
my $z1 = $ocean_thick;
my $vs1 = 0.5;
my $z2 = $z_max;
my $vs2 = 4.1;
my $z = 0.0;
my $ref_vmod_in = "../ref_vmod.in";

# Test model
my $vmod_out = "vmod.true";


# Make reference model 
open my $REF, ">", $ref_vmod_in or die "ERROR: cannot create $ref_vmod_in\n";
while ($z <= $z_max + 0.0001) {
    if ($z <= $ocean_thick) {
	print {$REF} "$z 1.5 0.01\n";
    }
    else {
	my $vs = ($vs2 - $vs1) / ($z2 - $z1) * ($z - $z1) + $vs1;
	my $vp = vs2vp($vs);
	print {$REF} "$z $vp $vs\n";
    }
    $z += $dz;
}
close $REF or die;

# Read reference model
my @vs_ref = @{read_ref_v($ref_vmod_in)};

# Make test model
open my $OUT, ">", $vmod_out or die "ERROR: Cannot create $vmod_out";
my $nlay = scalar(@dvs) + 1;
print {$OUT} "$nlay\n";
print {$OUT} "1.5 -1.0 1.0 $ocean_thick\n";
foreach my $i (0..$#dvs) {
    my $z_top 
	= $i == 0 ? $ocean_thick
	:           $zs[$i-1]
	;
    my $z_bot 
	= $i == $#dvs ? $z_max
	:               $zs[$i]
	;
    my $h = $z_bot - $z_top;
    my $vs = get_vs_ref($z_top, $z_bot, $dz, \@vs_ref);
    my $vp = vs2vp($vs);
    my $rho = vp2rho($vp);
    print {$OUT} "$vp $vs $rho $h\n";
}
close $OUT or die;

#-----------------------------------------------------------------------

sub read_ref_v {
    my $ref_vmod_in = $_[0];
    open my $REF,  "<", $ref_vmod_in or die;
    my (@vs_ref);
    while (my $line = <$REF>) {
	chomp $line;
	my @item = split q{ }, $line;
	push @vs_ref, $item[2];
    }
    close $REF or die;
    
    return \@vs_ref;
}

#-----------------------------------------------------------------------

sub get_vs_ref {
    my $z_top = $_[0];
    my $z_bot = $_[1];
    my $dz    = $_[2];
    my @vs_ref = @{$_[3]};

    my $zc = 0.5 * ($z_top + $z_bot);
    my $iz = int($zc / $dz + 0.5);
    my $n = @vs_ref;
    
    my $vs = $vs_ref[$iz];
    
    return $vs;
}

#-----------------------------------------------------------------------

sub vp2rho {
    my $vp = $_[0];
    my $rho = 1.6112 * $vp - 0.4721 * $vp**2 + 0.0671 * $vp**3
	- 0.0043 * $vp**4 + 0.000106 * $vp**5;
    return $rho;
}

#-----------------------------------------------------------------------

sub vs2vp {
    my $vs = $_[0];
    my $vp = 0.9409 + 2.0947 * $vs - 0.8206 * $vs**2
	+ 0.2683 * $vs**3 -  0.0251 * $vs**4;
    return $vp;
}
