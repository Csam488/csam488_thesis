+#!usr/bin/perl
#use warnings;
use strict;

#ranker script

package Whirley; # from http://www.perlmonks.org/?node_id=4943
sub new {
	my $type = shift;
        my $class = ref $type || $type;
        my $WHIRLEY_COUNT = -1;
        my @whirley = ("-", "\\", "|", "/");#map chr, qw/32 176 177 178 219 178 177 176/;
        my $self = sub {
		$WHIRLEY_COUNT = 0 if ++$WHIRLEY_COUNT == @whirley;
		return $whirley[$WHIRLEY_COUNT];
        };
        bless $self, $class;
}

#load spinner
my $whirley = new Whirley;
my $wirC = 1000;

#open stuff
my $name = "XX";
my $path = 'C:\Users\chris_000\Documents\Perl testing\\';
print ":::::::::::::::Starting:::::::::::::::\nRanking\n";
open TAR, $path."Tar.txt" or die "Can't open infile : $!"; 
my $Tar;
while (my $line = <TAR>){
	chomp $line;
	my @col = split /\t/, $line;
	$Tar = "$Tar:$col[0]";
}
$Tar = "$Tar:";
my $tot = 0;
my $a = 0;
my $b = 0;
my $c = 0;
my $d = 0;
my $e = 0;
my $f = 0;
my $g = 0;
open VAR, $path.$name.".vcf" or die "Can't open infile : $!";
open UTRT, ">".$path.$name."_UTR_TAR.vcf" or die "Can't open infile : $!"; #Note vestigial from older vestion
open CDT, ">".$path.$name."_COD_TAR.vcf" or die "Can't open infile : $!";
open CDD, ">".$path.$name."_COD_DAM.vcf" or die "Can't open infile : $!";
open CDB, ">".$path.$name."_COD_BRA.vcf" or die "Can't open infile : $!";
open CDO, ">".$path.$name."_COD_OTH.vcf" or die "Can't open infile : $!";
open PO, ">".$path.$name."_PO.vcf" or die "Can't open infile : $!";
while (my $line = <VAR>){
        if ($wirC == 1000){
		print STDERR "Working: ", $whirley->(), "\r";
		$wirC = 0;
	}else{
		$wirC++;
	}
	if ($line =~ /#/){
		print UTRT $line;
		print CDT $line;
		print CDD $line;
		print CDB $line;
		print CDO $line;
	}else{
		$tot++;
		my $sw = 0;
		my @col = split /\t/, $line;
		my $target = TARGET($col[7]);
		my $tissue = Tissue($col[7]);
		my $path = Path($col[7]);
		my $impact = Impact($col[7]);
		my $utr = UTR($col[7]);
		if ($utr == 1){
			if ($target == 1){
				print UTRT $line;
				$a++;
				$sw = 1;
			}
		}
		if ($impact > 0){
			if ($target == 1){
				print CDT $line;
				$c++;
				$sw = 1;
			}elsif ($tissue == 1 or $path == 1){
				print CDB $line;
				$e++;
				$sw = 1;
			}elsif ($impact == 2){
				print CDD $line;
				$d++;
				$sw = 1;
			}else{
				print CDO $line;
				$f++;
				$sw = 1;
			}
		}
		if ($sw == 0){
			#die $line;
		}
	}
}

print "\nUTR target $a\nCoding target $c\nCoding brain or path $e\nCoding dam $d\nCoding other $f\n\n\t>$tot\n";

sub Path{
	if ($_[0] =~ /hsa04080/ or $_[0] =~ /hsa04517/ or $_[0] =~ /hsa04724/ or $_[0] =~ /hsa04727/ or $_[0] =~ /hsa04725/ or $_[0] =~ /hsa04728/ or $_[0] =~ /hsa04726/ or $_[0] =~ /hsa04720/ or $_[0] =~ /hsa04730/ or $_[0] =~ /hsa04723/ or $_[0] =~ /hsa04721/ or $_[0] =~ /hsa04722/ or $_[0] =~ /hsa04040/){
		return 1;
	}else{
		return 0;
	}
}
sub TARGET{
	if ($_[0] eq "MAX"){
		return 2;
	}
	my @genes;
	if ($_[0] =~ /gSYM=([\w\d-]+)/){
		push @genes, $1;
	}
	if ($_[0] =~ /GL=([\w\d-]+)/){
		push @genes, $1;
	}
	foreach my $gene (@genes){
		if ($Tar =~ /:\Q$gene:/){
			return 1;
		}
	}
	return 0;
}
sub Tissue{
	my $r = 0;
	my $mat = 0;
	my $sp = 0;
	my $br = 0;
	if ($_[0] =~ /R_HPA=([\w \[\](),\/]*);/){
		my $temp = $1;
		$mat++;
		if ($temp =~ /en[rh]/){
			$sp++;
			if ($temp =~ /ereb[re]/ or $temp =~ /ippo/ or $temp =~ /ypotha/ or $temp =~ /ituit/ or $temp =~ /rain/ or $temp =~ /audate/){
				$br++;
			}
		}
	}if ($_[0] =~ /R_GTEx=([\w \[\](),\/]*);/){
		my $temp = $1;
		$mat++;
		if ($temp =~ /en[rh]/){
			$sp++;
			if ($temp =~ /ereb[re]/ or $temp =~ /ippo/ or $temp =~ /ypotha/ or $temp =~ /ituit/ or $temp =~ /rain/ or $temp =~ /audate/){
				$br++;
			}
		}
	}if ($_[0] =~ /R_FANTOM=([\w \[\](),\/]*);/){
		my $temp = $1;
		$mat++;
		if ($temp =~ /en[rh]/){
			$sp++;
			if ($temp =~ /ereb[re]/ or $temp =~ /ippo/ or $temp =~ /ypotha/ or $temp =~ /ituit/ or $temp =~ /rain/ or $temp =~ /audate/){
				$br++;
			}
		}
	}
	if ($mat == 0 or $br == 0){
		$r = 0;
	}elsif ($br >= ($mat / 2)){
		$r = 1;
	}elsif (($sp-$br) >= ($mat / 2)){
		$r = 0;
	}else{
		$r = 0;
	}
	return "$r";
}
sub Impact{
	if ($_[0] =~ /Cnsq=([\w:\-&_,\/\\]*);/){
		my $temp = $1;
		if ($temp =~ /stop/ or $temp =~ /start/ or $temp =~ /splice_[da]/ or $temp =~ /frame/ or $temp =~ /transcript_ablation/){
			return 2;
		}
	}
	if ($_[0] =~ /FG=([\w:_,\-\/\\]*);/){
		my $temp = $1;
		if ($temp =~ /stop/ or $temp =~ /splice-[da]/ or $temp =~ /frame/){
			return 2;
		}
	}
	if ($_[0] =~ /FD=([\w:_,\-\/\\]*);/){
		my $temp = $1;
		if ($temp =~ /stop/ or $temp =~ /splice-[da]/ or $temp =~ /frame/){
			return 2;
		}
	}
	if ($_[0] =~ /Cnsq=([\w:\-&_,\/\\]*);/){
		my $temp = $1;
		if ($temp =~ /miss/ or $temp =~ /protein/ or $temp =~ /splice/){
			return 1;
		}
	}
	if ($_[0] =~ /FG=([\w:_,\-\/\\]*);/){
		my $temp = $1;
		if ($temp =~ /miss/ or $temp =~ /protein/ or $temp =~ /splice/){
			return 1;
		}
	}
	if ($_[0] =~ /FD=([\w:_,\-\/\\]*);/){
		my $temp = $1;
		if ($temp =~ /miss/ or $temp =~ /protein/ or $temp =~ /splice/){
			return 1;
		}
	}
	return 0;
}
sub UTR{
	if ($_[0] =~ /Cnsq=([\w:&_\-,\/\\]*);/){
		my $temp = $1;
		if ($temp =~ /[Uu][Tt][Rr]/){
			return 1;
		}
	}
	if ($_[0] =~ /FG=([\w:_,\-\/\\]*);/){
		my $temp = $1;
		if ($temp =~ /[Uu][Tt][Rr]/){
			return 1;
		}
	}
	if ($_[0] =~ /FD=([\w:_,\-\/\\]*);/){
		my $temp = $1;
		if ($temp =~ /[Uu][Tt][Rr]/){
			return 1;
		}
	}
	return 0;
}
