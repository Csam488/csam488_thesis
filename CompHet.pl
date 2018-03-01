+#! usr/bin/perl
use strict;
use warnings;

################################################################
#NOTE: Old version of code as origonal was lost to file corruption 
################################################################
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
my $child = "XX";
my $mom = "XX";
my $dad = "XX";
my $thresh = "0.02";
my $c;
my $m;
my $d;
print ":::::::::::::::Starting:::::::::::::::\n";
open VAR, $path.$name.".vcf" or die "Can't open infile : $!";
my %cand;
my %chang;
use LWP::Simple;



while (my $line = <VAR>){
	if ($line =~ /^#/){
		if ($line =~ /#C/){
			chomp $line;
			my @col = split /\t/, $line;
			foreach my $x (9..$#col){
				if ($col[$x] eq $child){
					$c = $x;
				}				
				if ($col[$x] eq $mom){
					$m = $x;
				}				
				if ($col[$x] eq $dad){
					$d = $x;
				}
			}
		}	
	}else{
		chomp $line;
		my @col = split /\t/, $line;
		my $pos = "$col[0]:$col[1]";
		my @for = split /GT/, $col[8];
		my $GT = ($for[0] =~ s/:/:/g);
		my @forb = split /GQ/, $col[8];
		my $GQ = ($forb[0] =~ s/:/:/g);
		if (defined $GT and length $GT){
		}else{
			$GT = 0;
		}
		if (defined $GQ and length $GQ){
		}else{
			$GQ = 0;
		}
		my $cg;
		my $mg;
		my $dg;
		my $mq;
		my $cq;
		my $dq;
		foreach my $x (9..$#col){
			my @sam = split /:/, $col[$x];
			if ($x == $c){
				$cg = $sam[$GT];
				$cq = $sam[$GQ];
			}
			if ($x == $m){
				$mg = $sam[$GT];
				$mq = $sam[$GQ];
			}
			if ($x == $d){
				$dg = $sam[$GT];
				$dq = $sam[$GQ];
			}
		}
		if (defined $cq and length $cq and $cq ne "."){
		}else{
			$cq = 99;
		}if (defined $mq and length $mq and $mq ne "."){
		}else{
			$mq = 99;
		}if (defined $dq and length $dq and $dq ne "."){
		}else{
			$dq = 99;
		}
		my $PO;
		my @ca = split /\//, $cg;
		my @da = split /\//, $dg;
		my @ma = split /\//, $mg;
		my $fq = 1;
		my $no = 30;
		if ($ca[0] ne $ca[1] and $cq > $no and $mq > $no and $dq > $no){
			
			if ($mg eq "./." and $dg eq "./."){
			}elsif ($mg eq "./."){
				foreach my $ca (@ca){
					if ($ca != 0){
						if ($da[0] == $ca and $da[1] == $ca){
						}elsif ($da[0] == $ca or $da[1] == $ca){
							$PO = "F";
						}else{
							$PO = "Mp";
							$fq = freq ($col[7],$col[4],$ca,$col[3]);
						}
					}
				}
			}elsif ($dg eq "./."){
				foreach my $ca (@ca){
					if ($ca != 0){
						if ($ma[0] == $ca and $ma[1] == $ca){
						}elsif ($ma[0] == $ca or $ma[1] == $ca){
							$PO = "M";
							$fq = freq ($col[7],$col[4],$ca,$col[3]);
						}else{
							$PO = "Fp";
							$fq = freq ($col[7],$col[4],$ca,$col[3]);
						}
					}
				}
			}else{
				foreach my $ca (@ca){
					if ($ca != 0){
						if ($ma[0] == $ca and $ma[1] == $ca){
						}elsif ($da[0] != $ca and $da[1] != $ca and $ma[0] != $ca and $ma[1] != $ca){
							$PO = "De";
							my $temp = freq ($col[7],$col[4],$ca,$col[3]);
							if ($temp < $fq){
								$fq = $temp;
							}
						}elsif ($da[0] == $ca and $da[1] == $ca){
						}elsif ($ma[0] == $ca or $ma[1] == $ca){
							if ($da[0] == $ca or $da[1] == $ca){
							}else{
								if (defined $PO and length $PO){
									$PO = "Co";
									my $temp = freq ($col[7],$col[4],$ca,$col[3]);
									$fq = "$temp:$fq";
								}else{
									$PO = "M";
									$fq = freq ($col[7],$col[4],$ca,$col[3]);
								}
							}
						}elsif ($da[0] == $ca or $da[1] == $ca){
							if ($ma[0] == $ca or $ma[1] == $ca){
							}else{
								if (defined $PO and length $PO){
									$PO = "Co";
									my $temp = freq ($col[7],$col[4],$ca,$col[3]);
									$fq = "$fq:$temp";
								}else{
									$PO = "F";
									$fq = freq ($col[7],$col[4],$ca,$col[3]);
								}
							}
						}
					}
				}
			}
			if (defined $PO){
				$cg =~ /([123456789])/;
				my $thing = $1 - 1;
				my @alt = split //, $col[4];
				$chang{"$col[0]:$col[1]"} = "$col[3]-$alt[$thing]";
				my $cns = cons ($col[7]);
				my %genes;
				if($line !~ /GL=([\w\d-]*)/ and $line !~ /gSYM=([\w\d-]*)/){
				}else{
					if ($line =~ /GL=([\w\d-]*)/){
						$genes{$1} = 1;
					}
					if ($line =~ /gSYM=([\w\d-]*)/){
						$genes{$1} = 1;
					}
					my @genes = keys %genes;
					foreach my $gene (@genes){
						if (defined $cand{$gene}){
							$cand{$gene} = "$cand{$gene},$col[0]:$col[1];$PO;$fq;$cns";
						}else{
							$cand{$gene} = "$col[0]:$col[1];$PO;$fq;$cns";
						}
					}
				}
			}
		}
	}
	if ($wirC == 1000){
		print STDERR "Reading file: ", $whirley->(), "\r";
		$wirC = 0;
	}
	else{
		$wirC++;
	}
}
close VAR;
my $pass = 0;
my $sing = 0;
my $nocom = 0;
my $pr = 0;
print "\n:::::::::::::::Read Complete:::::::::::::::\n:::::::::::::::Calculating Compounds:::::::::::::::\n";
open TOT, ">".$path.$name."_COM_TOT.txt" or die "Can't open infile : $!";
open PASS, ">".$path.$name."_COM_PASS.txt" or die "Can't open infile : $!";
open PR, ">".$path.$name."_COM_PR.txt" or die "Can't open infile : $!";
open SD, ">".$path.$name."_COM_SD.txt" or die "Can't open infile : $!";
my %Exac;
open TAR, $path."Tar.txt" or die "Can't open infile : $!"; 
my $Tar;
while (my $line = <TAR>){
	chomp $line;
	my @col = split /\t/, $line;
	$Tar = "$Tar:$col[0]";
}
close TAR;
$Tar = "$Tar:";
my @genes = keys %cand;
my $excnt = 0;
foreach my $gene (@genes){
	if ($Tar =~ /\Q$gene:/){
	my @entries = split /,/, $cand{$gene};
	if ($#entries){
		my $M =0;
		my $F =0;
		my $B =0;
		foreach my $entry (@entries){
			if ($entry =~ /;M/){
				$M++;
			}elsif ($entry =~ /;F/){
				$F++;
			}elsif ($entry =~ /;[CD]/){
				$B++;
			}else{
				die "Something is bad with origin calculation";
			}
		}
		my $COMBO = ($F*$M)+($B*($M+$F+$B));
		if ($COMBO and $gene !~ /^HLA/ and $gene !~ /OR[0-9]/){
			my @bit;
			foreach my $x (0..($COMBO-1)){
				$bit[$x] = 0;
			}
			my $flag = join "", @bit;
			foreach my $x (0..$#entries){
				$entries[$x] = "$entries[$x];$flag";
			}
			my $b = 0;
			my @filt;
			my @HH;
			my $f = 0;
			foreach my $x (0..($#entries-1)){
				foreach my $y (($x+1)..$#entries){
					if (($entries[$x] =~ /;[FCD]/ and $entries[$y] =~ /;[MCD]/) or ($entries[$y] =~ /;[FCD]/ and $entries[$x] =~ /;[MCD]/)){
						my @first = split /;/, $entries[$x];
						my @fbit = split //, $first[$#first];
						my @second = split /;/, $entries[$y];
						my @sbit = split //, $second[$#second];
						$fbit[$b] = 1;
						$sbit[$b] = 1;
						if ($b >= $COMBO){
							die "ERROR\n@entries\n$F\n$M\n$B\n$COMBO\t$b\n"
						}
						$b++;
						my $fbit = join "", @fbit;
						my $sbit = join "", @sbit;
						$first[$#first] = $fbit;
						$second[$#second] = $sbit;
						my $f = join ";", @first;
						my $s = join ";", @second;
						$entries[$x] = $f;
						$entries[$y] = $s;
						if ($first[2] =~ /:/ and $second[2] =~ /:/){
							 my @fthings = split /:/, $first[2];
							 my @fs = sort @fthings;
							 my @sthings = split /:/, $second[2];
							 my @ss = sort @sthings;
							 if (($fs[1]*$ss[1]) <= $thresh){
							 	push @filt, "$x:$y";
							 	if ($first[3] =~ /^h/ and $second[3] =~ /^h/){
							 		push @HH, "$x:$y";
							 	}
							 }
						}elsif ($first[2] =~ /:/){
						        my @things = split /:/, $first[2];
						        if ($second[1] =~ /M/){
						        	if (($things[1]*$second[2]) <= $thresh){
						        		push @filt, "$x:$y";
						        		if ($first[3] =~ /^h/ and $second[3] =~ /^h/){
							 		push @HH, "$x:$y";
							 	}
						        	}
						        }elsif ($second[1] =~ /F/){
						        	if (($things[0]*$second[2]) <= $thresh){
						        		push @filt, "$x:$y";
						        		if ($first[3] =~ /^h/ and $second[3] =~ /^h/){
							 		push @HH, "$x:$y";
							 	}
						        	}
						        }else{
						        	my @sthings = sort @things;
						        	if (($sthings[0]*$second[2]) <= $thresh){
						        		push @filt, "$x:$y";
						        		if ($first[3] =~ /^h/ and $second[3] =~ /^h/){
							 		push @HH, "$x:$y";
							 	}
						        	}
						        }
						}elsif ($second[2] =~ /:/){
							my @things = split /:/, $second[2];
						        if ($first[1] =~ /M/){
						        	if (($things[1]*$first[2]) <= $thresh){
						        		push @filt, "$x:$y";
						        		if ($first[3] =~ /^h/ and $second[3] =~ /^h/){
							 		push @HH, "$x:$y";
							 	}
						        	}
						        }elsif ($first[1] =~ /F/){
						        	if (($things[0]*$first[2]) <= $thresh){
						        		push @filt, "$x:$y";
						        		if ($first[3] =~ /^h/ and $second[3] =~ /^h/){
							 		push @HH, "$x:$y";
							 	}
						        	}
						        }else{
						        	my @sthings = sort @things;;
						        	if (($sthings[0]*$first[2]) <= $thresh){
						        		push @filt, "$x:$y";
						        		if ($first[3] =~ /^h/ and $second[3] =~ /^h/){
							 		push @HH, "$x:$y";
							 	}
						        	}
						        }
						}elsif (($first[2]*$second[2]) <= $thresh){
							my $ExAC = 0;
							#die $ExAC;
							if ($first[3] =~ /^[hm][ie]/ and $second[3] =~ /^[hm][ie]/){
								if ($first[2]*$first[2] > $thresh){
									$ExAC = ExAC($first[0],$first[1]);
								}elsif ($second[2]*$second[2] > $thresh){
									$ExAC = ExAC($second[0],$second[1]);
								}
								if ($ExAC <= $thresh){
									push @filt, "$x:$y";
									if ($first[3] =~ /high/ and $second[3] =~ /high/){
										push @HH, "$x;$y";
									}
								}
							}
						}
					}
					if ($wirC == 1000){
				print STDERR "Working: ", $whirley->(), "\r";
				$wirC = 0;
			}
			else{
				$wirC++;
			}
				}
			}
			
			if (defined $filt[0]){
				my %filt;
				foreach my $x (0..$#filt){
					my @thing = split /:/, $filt[$x];
					foreach my $thing (@thing){
						if (defined $filt{$thing}){
							my @frbit = split //, $filt{$thing};
							$frbit[$x] = 1;
							my $frbit = join "", @frbit;
							$filt{$thing} = $frbit;
						}else{
							my @frbit;
							foreach my $y (0..$#filt){
								push @frbit, "0";
							}
							$frbit[$x] = 1;
							my $frbit = join '', @frbit;
							$filt{$thing} = $frbit;
						}
					}
				}
				my @f = keys %filt;
				foreach my $f (@f){
					$entries[$f] = "$entries[$f];$filt{$f}";
				}
			}
			foreach my $x (0..$#entries){
				print TOT "$gene\t$entries[$x]\n";
				my @thing = split /;/, $entries[$x];
				if (defined $thing[5]){
					print PASS "$gene\t$entries[$x]\n";
					$pass++;
					my $HH = join ";", @HH;
					$HH = ";$HH";
					if ($HH =~ /;\Q$x/){
						print PR "$gene\t$entries[$x]\n";
						$pr++;
					}else{
						print SD "$gene\t$entries[$x]\n";
					}					
				}
			}
		}else{
			$nocom++;
			my @bit;
			foreach my $x (0..($COMBO-1)){
				$bit[$x] = 0;
			}
			my $flag = join "", @bit;
			foreach my $entry (@entries){
				print TOT "$entry;$flag\n";
			}
		}	
	}else{
		$sing++;
		print TOT "$gene\t$cand{$gene}\n";
	}
}
}
close TOT;
close PASS;
close PR;
close SD;
print "\n:::::::::::::::Calculations Complete:::::::::::::::\n:::::::::::::::Generating OutFiles:::::::::::::::\n";
my @files = qw ( PR SD );
foreach my $file (@files){
	open FIL, $path.$name."_COM_".$file.".txt" or die "Can't open infile : $!";
	open OUT, ">".$path.$name."_COM_".$file.".vcf" or die "Can't open infile : $!";
	my @chrom;
	my @pos;
	my @PO;
	my @MbF;
	my @FbF;
	while (my $line = <FIL>){
		chomp $line;
		my @ent = split /\t/, $line;
		my @entb = split /;/, $ent[1];
		my @entc = split /:/, $entb[0];
		push @chrom, $entc[0];
		push @pos, "$entc[0]\t$entc[1]";
		push @PO, $entb[1];
		push @MbF, $entb[4];
		push @FbF, $entb[5];
	}
	if ($#chrom =! $#pos){
		die "ERROR::\tBad split\n";
	}
	close FIL;
	open VAR, $path.$name.".vcf" or die "Can't open infile : $!";
	while (my $line = <VAR>){
		if ($wirC == 1000){
				print STDERR "Working: ", $whirley->(), "\r";
				$wirC = 0;
			}
			else{ 
				$wirC++;
			}
		if ($line =~ /^##/){
			print OUT $line;
		}elsif ($line =~ /^#/){
			print OUT '##INFO=<ID=PO,Number=1,Type=String,Description="Parental origin">',"\n";
			print OUT '##INFO=<ID=MbF,Number=1,Type=String,Description="Bit flag of matches per gene">',"\n";
			print OUT '##INFO=<ID=Fbf,Number=1,Type=String,Description="Bit flag of matches passing filter per gene">',"\n";
			print OUT $line;
		}else{
			foreach my $x (0..$#pos){
				if ($line =~ /\Q$pos[$x]/){
					my @col = split /\t/, $line;
					$col[7] = "$col[7];PO=$PO[$x];MbF=$MbF[$x];FbF=$FbF[$x]";
					my $out = join "\t", @col;
					print OUT $out;
				}
			}
		}
	}
}
					
				
		 
print "\n:::::::::::::::Complete:::::::::::::::\n\t";
print $#genes +1;
print " Genes analysed\n\t$sing genes had one candidate\n\t$nocom genes had no combinations\n\n\t$pass entries passed the filter\n\t$pr prime candidates\n";





			
sub ExAC{
	my @ref = split /:/, $_[0];
	my $ref = "$ref[0]-$ref[1]";
	if (defined $Exac{"$_[0]"}){
		#return $Exac{"$_[0]"};
	}
	#print "$ref\n";
	my $url = 'http://exac.broadinstitute.org/variant/'."$ref-$chang{$_[0]}";
	#$escnt++;
	#print "$ecnt","\r";
	sleep 30;
	my $content = get $url;
	my $ho = 0;
	if (defined $content){
		my $table = "NAN";
		if($content =~ /(<table id="frequency_table"[\w\W]*?<\/table>)/){
			$table = $1;
		}
		my @entries;
		while ($table =~ s/<tr>([\w\W]*?)<\/tr>//){
			push @entries, $1;
		}
		foreach my $n (1..($#entries-1)){
			my @values;
			while ($entries[$n] =~ s/>([\d.]+)<//){
				push @values, $1;
				#print "$1\n"
			}
			#print "=================\n";
			if ($#values == 3){
				my $E_freq = $values[2]/($values[1]/2);
				if ($E_freq > $ho){
					$ho = $E_freq;
				}
			}
		}
	}
	if ($ho > 0.01){
		#die "$ref-$chang{$_[0]}\n$ho\n";
		#print "$ref-$chang{$_[0]}\n$ho\n";
	}
	$Exac{$_[0]} = $ho;
	return $ho;	
	
}	


sub cons{
	if ($_[0] =~ /stop[-_]/ or $_[0] =~ /frameshift/ or $_[0] =~ /splice-[da]/ or $_[0] =~ /start[-_][lL]/){
		return "high";
	}elsif ($_[0] =~ /missense/ or $_[0] =~ /inframe/){
		return "med";
	}elsif ($_[0] =~ /UTR/ or $_[0] =~ /utr/){
		return "utr";
	}else{
		return "low";
	}
}



sub freq{
	my @bases = qw( FQ_GEXAC FQ_GESP FQ_1KG );
	my @info = split /;/, $_[0];
	my @alt = split /,/, $_[1];
	my $fpos = $_[2]-1;
	my $base = $alt[$fpos];
	my $temp = 0;
	if (length $base > length $_[3]){
		$base =~ s/^\Q$_[3]//;
		$base = "ins$base";
	}elsif (length $base < length $_[3]){
		my $t = $_[3];
		$t =~ s/^\Q$base//;
		$base = "del$t";
	}elsif ($base eq "*"){
		$base = "del$_[3]";
	}
	foreach my $i (@info){
		foreach my $b (@bases){
			if ($i =~ s/\Q$b=//){
				if (length $i){
					my @freq = split /\//, $i;
					foreach my $f (@freq){
						if ($f =~ /^\Q$base:/){
							$f =~ s/\Q$base://;
							if ($temp > $f or $temp == 0){
								$temp = $f;
							}
						}
					}
				}
			}
		}
	}
	return "$temp";
	
}
	
