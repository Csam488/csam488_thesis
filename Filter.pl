+#! usr/bin/perl
use strict;
#use warnings;


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
close FREQ;
my @bases = qw( FQ_GEXAC FQ_1KG FQ_GESP );
my $name = "XX";
my $path = 'C:\Users\chris_000\Documents\Perl testing\\';
my $child = "XX";
my $freq = 0.02;
my @names = qw( OR HLA );
my %renam;

open TAR, $path."Tar.txt" or die "Can't open infile : $!"; 
my $Tar = ':';
while (my $line = <TAR>){
	chomp $line;
	my @col = split /\t/, $line;
	$Tar = "$Tar:$col[0]";
}
close TAR;
$Tar = "$Tar:";
my @step;
my @rem;
foreach my $x (0..9){
	$rem[$x] = 0;
}
my %gen;
open GEN, $path."Human_genes.txt" or die "Can't open genes : $!";
while (my $line = <GEN>){
	my @col = split /\t/, $line;
	$gen{$col[2]} = $col[9];
}
close GEN;

my %freq;
open FREQ, $path."Freq.txt" or die "Can't open freq : $!"; 
while (my $line = <FREQ>){
	chomp $line;
	my @col = split /\t/, $line;
	$freq{$col[0]} = $col[1];
}
print ":::::::::::::::Starting:::::::::::::::\nGene Black List:\n";
foreach my $name (@names){
	print "\t-$name\n";
}
print "Filtering frequency:\n\t$freq\n";
open VAR, $path.$name.".vcf" or die "Can't open infile : $!";
open OUT, ">".$path.$name."_Filt.vcf" or die "Can't open outfile : $!";
open OUB, ">".$path.$name."_PDel.vcf" or die "Can't open outfile : $!";
my $ch;
my $ct = 0;
my $rm = 0;
my $non = 0;
while (my $line = <VAR>){
	if ($line =~ /^#/){
		if ($line =~ /^#CHROM/){
			chomp $line;
			my @col = split /\t/,$line;
			my( $i )= grep { $col[$_] eq $child } 0..$#col;
			$ch = $i;
			$line = "$line\n";
			print OUT '##INFO=<ID=',"FILT",',Number=.,Type=String,Description="Alleles that did not pass frequency filter">',"\n";
		}
		print OUT $line;
		print OUB $line;
	}else{
		my %NO;
		my $freqtst = 1;
		my $switch = 1;
		while ($switch == 1){
			my $ent = $line;
			chomp $ent;
			my @col = split /\t/, $ent;
			my @format = split /:/, $col[8];
			my @child = split /:/, $col[$ch];
			my( $g )= grep { $format[$_] eq "GT" } 0..$#format;
			my( $q )= grep { $format[$_] eq "GQ" } 0..$#format;
			if ($child[$g] eq './.' or $child[$g] eq '0/0'){
				if ($child[$g] eq './.' ){
					print OUB $line;
				}
				$freqtst = 0;
				$switch = 0;
				$step[0]++;
			}
			last if ($switch == 0);
			$rem[0]++;
			if ($child[$q] < 30){
				$freqtst = 0;
				$switch = 0;
			}
			last if ($switch == 0);
			$rem[1]++;
			foreach my $base (@bases){
				if ($line =~ /\Q$base=/){
					my $result = freq($base,$line);
					if ($result eq "K"){
						$freqtst = 0;
						$switch = 0;	
					}elsif ($result ne ">" and length $result and defined $result){
						$NO{$result} = 1;
					}
				}
				last if ($switch == 0);
			}
			last if ($switch == 0);
			my @NO = keys %NO;
			push @NO,"0";
			if ($#NO > 0){
				my $result = comb($line,@NO);
				if ($result eq "1"){
					$freqtst = 0;
					$switch = 0;	
				}
			}
			last if ($switch == 0);
			$rem[2]++;
			if ($line =~ /Cnsq=/ or $line =~ /FG=/ or $line =~ /FD=/){
				my $result = type($line);
				if ($result == 1){
					$freqtst = 0;
					$switch = 0;	
				}
			}
			last if ($switch == 0);
			$rem[3]++;
			if ($line !~ /Cnsq=/ and $line !~ /FG=/ and $line !~ /FD=/){
				$freqtst = 0;
				$switch = 0;
			}
			last if ($switch == 0);
			$rem[4]++;
			if ($line =~ /;IM=[PO]/ or $line =~ /;IM=DE_0/){
				$freqtst = 0;
				$switch = 0;
			}
			last if ($switch == 0);
			$rem[5]++;
			foreach my $nam(@names){
				if ($switch == 0){
				}else{
					if ($line =~ /GL=\Q$nam/ and $line =~ /gSYM=\Q$nam/ and $line !~ /IM=O/){
						$line =~ /GL=([\w\d-]*)/;
						last if ($1 =~ /OR[A-Z]/);
						$renam{$1}=1;
						$line =~ /gSYM=([\w\d-]*)/;
						last if ($1 =~ /OR[A-Z]/);
						$renam{$1}=1;
						$freqtst = 0;
						$switch = 0;
					}
					elsif ($line !~ /GL=/ and $line =~ /gSYM=\Q$nam/ and $line !~ /IM=O/){
						$line =~ /gSYM=([\w\d-]*)/;
						last if ($1 =~ /OR[A-Z]/);
						$renam{$1}=1;
						$freqtst = 0;
						$switch = 0;
					}
					elsif ($line =~ /GL=\Q$nam/ and $line !~ /gSYM=/ and $line !~ /IM=O/){
						$line =~ /GL=([\w\d-]*)/;
						last if ($1 =~ /OR[A-Z]/);
						$renam{$1}=1;
						$freqtst = 0;
						$switch = 0;
					}
				}
			}
			last if ($switch == 0);
			$rem[6]++;
			my @genes;
			if ($line =~ /gSYM=([\w\d\.-]+)/){
				push @genes, $1;
				}
			if ($line =~ /GL=([\w\d\.-]+)/){
					push @genes, $1;
			}
			foreach my $gene (@genes){
				if (defined $gen{$gene}){
					if ($gen{$gene} eq "protein-coding"){
						if ($gene =~ /\./){
							#die $gene;
						}
					}else{
						$freqtst = 0;
						$switch = 0;
					}
				}elsif ($gene !~ /\./){
				}else{
					$freqtst = 0;
					$switch = 0;
				}
			}
			last if ($switch == 0);
			$rem[7]++;
			if ($#genes == 1){
				if ($genes[0] eq $genes[1]){
					if ($line =~ /IMP=[ATGC]*:MODIFIER/ or $line =~ /IMP=[ATGC]*:LOW/){
						if ($line =~ /near-splice/ or $line =~ /splice_region_variant/){
							$line =~ /DSP=([\d]*)/;
							if ($1 <= 99){
							}else{
							$freqtst = 0;
							$switch = 0;
							}
						}else{
							$freqtst = 0;
							$switch = 0;
						}
					}
				}
			}
			last if ($switch == 0);
			$rem[8]++;
			my @col = split /\t/, $line;
			my @alt = split /,/, $col[4];
		
		
			if ($col[$ch] =~ /^0\/1/){
				if (defined $freq{"$col[0]-$col[1]-$col[3]-$alt[0]"}){
					if ($freq{"$col[0]-$col[1]-$col[3]-$alt[0]"} ne "pass"){
						$freqtst = 0;
						$switch = 0;
					}
				}
			}elsif ($col[$ch] =~ /^1\/1/){
				if (defined $freq{"$col[0]-$col[1]-$col[3]-$alt[0]"}){
					if ($freq{"$col[0]-$col[1]-$col[3]-$alt[0]"} ne "pass" and $freq{"$col[0]-$col[1]-$col[3]-$alt[0]"} ne "hom"){
						$freqtst = 0;
						$switch = 0;
					}
				}
			}else{
				#print "other\n";
			}
			last if ($switch == 0);
			$rem[9]++;
			
			
			$switch = 0;			
		}
		if ($freqtst){
			my @NO = keys %NO;
			if(length $NO[0] and defined $NO[0]){
				chomp $line;
				my @col = split /\t/, $line;
				foreach my $x (0..6){
					print OUT "$col[$x]\t";
				}
				print OUT "$col[7];FILT=";
				if ($#NO == 0){
					print OUT "$NO[0]";
				}else{
					my $y = $#NO-1;
					foreach my $x (0..$y){
						print OUT $NO[$x],'/';
					}
					print OUT $NO[$#NO];
				}
				foreach my $x (8..$#col){
					print OUT "\t$col[$x]";
				}
				print OUT "\n";
			}else{
				print OUT $line;
			}
		}else{
			$rm++;
		}
		if ($wirC == 1000){
			print STDERR "Working: ", $whirley->(), "\r";
			$wirC = 0;
		}
		else{
			$wirC++;
		}
		$ct++;
	}
}
my $remain = $ct-$rm;
my $devhom = $remain-$non;
close OUB;
print "\n:::::::::::::::Completed:::::::::::::::\n\tRemoved $rm of $ct entries\n\t$remain remaining\nSteps\n";
my @nam = qw( Not_0/0_./. Qual Freq Consequnetial Annotated_to_Gene Mode Not_in_excluded_genes Prot_cod Unused_TranScpt In_House );
my $tot = 0;
foreach my $x(0..$#rem){
	print "$nam[$x]\t$rem[$x]\n";
}
my @prem = keys %renam;
my @remo = sort { $a <=> $b || $a cmp $b} @prem;
print "Genes removed:\n";
foreach my $remo (@remo){
	print "\t$remo\n";
}
print "\n\n";


sub type{
	my @col = split /\t/, $_[0];
	my $r = 1;
	if ($col[7] =~ /stop[-_]gain/ or $col[7] =~ /start[-_]loss/ or $col[7] =~ /stop[-_]loss/ or $col[7] =~ /splice[-_][da]/ or $col[7] =~ /frameshift/ or $col[7] =~ /missense/ or $col[7] =~ /inframe_/){
		$r = 0;
	}
	if ($r == 1){
		#print "\t@cons\n\n";
	}else {
		#print "@cons\n\n";
	}
	return $r;
	
}

sub comb{
	my @col = split /\t/, $_[0];
	my @format = split /:/, $col[8];
	my ( $i )= grep { $format[$_] eq "GT" } 0..$#format;
	my @cld = split /:/, $col[$ch];
	my @gt = split /\//, $cld[$i];
	my @alt = split /,/, $col[4];
	unshift @alt, $col[3];
	my $first;
	my $second;
	foreach my $x (1..$#_){
		if ($gt[0] eq $_[$x] or $gt[0] eq $col[3]){
			$first = ">";
		}
		if ($gt[1] eq $_[$x]  or $gt[0] eq $col[3]){
			$second = ">";
		}
	}
	if ($first eq $second and $second eq ">"){
		return 1;
	}else{
		return 0;
	}
}
	
sub freq{
	my @col = split /\t/, $_[1];
	my @format = split /:/, $col[8];
	my ( $i )= grep { $format[$_] eq "GT" } 0..$#format;
	my @cld = split /:/, $col[$ch];
	my @gt = split /\//, $cld[$i];
	if ($gt[0] eq "." or $gt[1] eq "."){
		return "K"
	} 
	my @info = split /;/, $col[7];
	my ( $j )= grep { $info[$_] =~ /^\Q$_[0]=/ } 0..$#info;
	$info[$j] =~ s/\Q$_[0]=//;
	my @FQ = split /\//, $info[$j]; 
	my @alt = split /,/, $col[4];
	unshift @alt, $col[3];
	my $first = @alt[$gt[0]];
	my $second = @alt[$gt[1]];
	my $fl = length($first);
	my $sl = length($second);
	my $rl = length($col[3]);
	my $f = 0;
	my $s = 0;
	if ($fl == $rl){
		if ($info[$j] !~ /\Q$first:/){
			#return "0";
		}else{
			my( $k )= grep { $FQ[$_] =~ /^\Q$first:/ } 0..$#FQ;
			$FQ[$k] =~ s/^\Q$first://;
			if ($FQ[$k] > $freq){
				$f = 1;
			}
		}
	}else{
		if ($fl > $rl){
			my $in = $first;
			$in =~ s/^\Q$col[3]//;
			my $in = "ins$in";
			if ($info[$j] !~ /\Q$in:/){
				#return "0";
			}else{
				my( $k )= grep { $FQ[$_] =~ /^\Q$in:/ } 0..$#FQ;
				$FQ[$k] =~ s/^\Q$in://;
				if ($FQ[$k] > $freq){
					$f = 1;
				}
			}
		}else{
			my $del = $col[3];
			$del =~ s/^\Q$first//;
			my $in = "del$del";
			if ($info[$j] !~ /\Q$in:/){
				#return "0";
			}else{
				my( $k )= grep { $FQ[$_] =~ /^\Q$in:/ } 0..$#FQ;
				$FQ[$k] =~ s/^\Q$in://;
				if ($FQ[$k] > $freq){
					$f = 1;
				}
			}
		}
	}
			
	if ($first eq $second){
		$s = $f;
	}else{
		if ($sl == $rl){
			if ($info[$j] !~ /\Q$second:/){
				#return "0";
			}else{
				my( $l )= grep { $FQ[$_] =~ /^\Q$second:/ } 0..$#FQ;
				$FQ[$l] =~ s/^\Q$second://;
				if ($FQ[$l] > $freq){
					$s = 1;
				}
			}
		}else{
			if ($sl > $rl){
				my $in = $second;
				$in =~ s/^\Q$col[3]//;
				my $in = "ins$in";
				if ($info[$j] !~ /\Q$in:/){
					#return "0";
				}else{
					my( $k )= grep { $FQ[$_] =~ /^\Q$in:/ } 0..$#FQ;
					$FQ[$k] =~ s/^\Q$in://;
					if ($FQ[$k] > $freq){
						$s = 1;
					}
				}
			}else{
				my $del = $col[3];
				$del =~ s/^\Q$second//;
				my $in = "del$del";
				if ($info[$j] !~ /\Q$in:/){
					#return "0";
				}else{
					my( $k )= grep { $FQ[$_] =~ /^\Q$in:/ } 0..$#FQ;
					$FQ[$k] =~ s/^\Q$in://;
					if ($FQ[$k] > $freq){
						$s = 1;
					}
				}
			}
		}
	}
	if ($s == 1 and $f == 1){
		return "K";
	}else{
		if ($s == 1 or $f == 1){
			if ($s == 1){
				my ( $r )= grep { $alt[$_] eq $second } 0..$#alt;
				return "$r";
			}else{
				my ( $r )= grep { $alt[$_] eq $first } 0..$#alt;
				return "$r";
			}
		}else{
			return ">";
		}
	}
}
