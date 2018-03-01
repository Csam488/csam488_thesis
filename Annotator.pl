#! usr/bin/perl
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

#Unified annotator program
#For caculating allele frequencies, giving trios modes of inheritance

my $manual = 0; #set manual entry or not
my $screen = 1; #set screen output

my @bases;
my $name;
my $path;
my $child;
my $mother;
my $father;
my $male;
my $freqsw;
my $modesw;
my $ext;
my $cto = 0;
if ($manual < 1){
	@bases = qw(GEXAC GESP); #array of id's of allele counts to be converted to frequencies
	$name = "XXX"; #file name
	$ext = ""; #fire extention
	$path = 'C:\Users\chris_000\Documents\Perl testing\\'; #file path
	$child = "XX"; #child code
	$mother = "XX"; #mother code
	$father = "XX"; #father code
	$male = 1; #male state of child (0 = female, 1 = male)
	$modesw = 1;
	$freqsw = 1;
}else{ #take manual input
	print "What is file name: ";
	$name = <>;
	print "What is file extention: ";
	$ext = <>;
	print "\nWhat is file path: ";
	$path = <>;
	print "\nWhat is the childs identifier: ";
	$child = <>;
	print "\nIs the child male (1=yes/0=no): ";
	$male = <>;
	print "\nWhat is the mothers identifier: ";
	$mother = <>;
	print "\nWhat is the fathers identifier: ";
	$father = <>;
	print "\nAre allele frequencies to be calculated (1=yse/0=no): ";
	$freqsw = <>;
	if ($freqsw){
		print "\nWhat are the ID's of the allele counts to be converted to frequencies (separated by spaces): ";
		my $temp = <>;
		@bases = split / /, $temp;
	}
	print "\nAre modes of inheitance to be predicted (1=yes/0=no): ";
	$modesw = <>;
}
my $whirley = new Whirley;
my $wirC = 1000;


#read file and dtermine number of entries	
print ":::::::::::::::Starting:::::::::::::::\n\tOpening infile\n";
my $entcnt = 0;
my @file;
my %PAT;
open VAR, $path.$name.$ext.".vcf" or die "Can't open infile : $!";
print "\tReading number of entries\n";
while (my $l = <VAR>){
	if ($l !~ /^#/){
		$entcnt++;
		if ($l =~ /gSYM=([\w\d-]*)/){
			$PAT{$1}="NA\tNA\tNA\tNA\tNA";
		}if ($l =~ /GL=([\w\d-]*)/){
			$PAT{$1}="NA\tNA\tNA\tNA\tNA";
		}
	}
	if ($screen){
		if ($wirC == 1000){
			print STDERR "Reading file: ", $whirley->(), "\r";
			$wirC = 0;
		}
		else{
			$wirC++;
		}
	}
	push @file,$l;
}
close VAR;
print "\n\t$entcnt entries found\n:::::::::::::::Retrieving HPA Info:::::::::::::::\n";
my @genes = keys %PAT;
my $PAT_ct = 0;
my %pos;
open IND, 'C:\Users\chris_000\Documents\Perl testing\Database.indx' or die ("Can't open infile : $!");
while (my $in =<IND>){
	chomp $in;
	my @ind = split /\t/, $in;
	$pos{$ind[0]} = $ind[1];
}
close IND;
use Tie::File;
tie my @IDX, 'Tie::File', 'C:\Users\chris_000\Documents\Perl testing\Database.tdt' or die ("Can't open infile : $!"); #######Expression_database
foreach my $gene (@genes){
	if ($screen){
		if ($wirC >= 1000){
			print STDERR "Reading entry: ", $whirley->(), "\r";
			$wirC = 0;
		}
		else{
			$wirC++;
		}
	}
	my $PAT_ent = $IDX[$pos{$gene}];
	if ($PAT_ent =~ /^\Q$gene/){
		chomp $PAT_ent;
		my @PAT_col = split /\t/, $PAT_ent;
		my @PAT_pcl;
		if ($PAT_col[13] =~ /Enzyme/){
			push @PAT_pcl, "Enzyme";
		}
		if ($PAT_col[13] =~ /CD/){
			push @PAT_pcl, "CD_marker";
		}
		if ($PAT_col[13] =~ /Blood/){
					push @PAT_pcl, "Blood_grp_antigen";
		}
		if ($PAT_col[13] =~ /Nuclear/){
			push @PAT_pcl, "Nuclear_receptor";
		}
		if ($PAT_col[13] =~ /Transport/){
			push @PAT_pcl, "Transporter";
		}
		if ($PAT_col[13] =~ /Ribosomal/){
			push @PAT_pcl, "Ribosomal_protein";
		}
		if ($PAT_col[13] =~ /coupled/){
			push @PAT_pcl, "G_protein_cpld_receptor";
		}
		if ($PAT_col[13] =~ /Voltage/){
			push @PAT_pcl, "Voltage_gated_ion_chnl";
		}
		if ($PAT_col[13] =~ /memberane/){
			push @PAT_pcl, "Prd_membrane_protein";
		}
		if ($PAT_col[13] =~ /secreted/){
			push @PAT_pcl, "Prd_secreted_protein";
		}
		if ($PAT_col[13] =~ /intracel/){
			push @PAT_pcl, "Prd_intracellular_protein";
		}
		if ($PAT_col[13] =~ /Plasma/){
			push @PAT_pcl, "Plasma_protein";
		}
		if ($PAT_col[13] =~ /Transcri/){
			push @PAT_pcl, "Transcription_factor";
		}
		if ($PAT_col[13] =~ /polym/){
			push @PAT_pcl, "RNA_pol_related";
		}
		if ($PAT_col[13] =~ /RAS/){
			push @PAT_pcl, "RAS_pth";
		}
		if ($PAT_col[13] =~ /Citric/){
			push @PAT_pcl, "Citric_acid_cycle_related";
		}
		my $PAT_class = join ",", @PAT_pcl;
		my @PAT_psum = split /;/, $PAT_col[6];
		my @PAT_br;
		foreach my $pat (@PAT_psum){
			if ($pat =~ /^cerebral cortex:([\d.]*)/){
				push @PAT_br, "Cerebral_cortex:$1";
			}elsif ($pat =~ /^hippocampus:([\d.]*)/){
				push @PAT_br, "Hippocampus:$1";
			}elsif ($pat =~ /^caudate:([\d.]*)/){
				push @PAT_br, "Caudate:$1";
			}elsif ($pat =~ /^cerebellum:([\d.]*)/){
				push @PAT_br, "Cerebellum:$1";
			}
		}
		my $PAT_pbr = join ",", @PAT_br;
		my @PAT_sum = split ";", $PAT_col[9];
		my @P_ta = split "=", $PAT_sum[0];
		my $PAT_hpa = $P_ta[1];
		my @P_tb = split "=", $PAT_sum[1];
		my $PAT_gtex = $P_tb[1];
		my @P_tc = split "=", $PAT_sum[2];
		my $PAT_fan = $P_tc[1];
		$PAT{$gene} = "$PAT_class\t$PAT_pbr\t$PAT_hpa\t$PAT_gtex\t$PAT_fan";
	}
}

close PAT;



print "\n\tComplete\n:::::::::::::::Annotating:::::::::::::::\n";
#report what annorating on
if ($freqsw){
	print "\tFrequency\n";
}
if ($modesw){
	print "\tMode\n";
}
print "\tHPA data\n";

#open files
open VAR, $path.$name.$ext.".vcf" or die "Can't open infile : $!";
open OUT, ">".$path.$name."_ANN.vcf" or die "Can't open outfile : $!";

my $PR;
my $MO;
my $FR;
my $cnt = 0;
my $cntD = 0;
my $cntR = 0;
my $cntX = 0;
my $cntE = 0;
my $cntP = 0;
my $cntL = 0;
my $cntO = 0;
my $freqCt = 0;
foreach my $line (@file){
	if ($line =~ /^##/){
		print OUT $line;
	}elsif ($line =~ /^#C/){
		#print info headders
		if ($freqsw){
			foreach my $base(@bases){
				print OUT '##INFO=<ID=',"FQ_$base",',Number=.,Type=String,Description="Allele Frequencies based on ',"$base",'">',"\n";
			}
		}
		if ($modesw){
			print OUT '##INFO=<ID=IM,Number=1,Type=String,Description="Inheritance mode based on genotypes of trio">',"\n";
			print OUT '##INFO=<ID=PCl,Number=1,Type=String,Description="Protein class HPA derived">',"\n";
			print OUT '##INFO=<ID=Pbr,Number=1,Type=String,Description="Protein expression in brian, HPA derived">',"\n";
			print OUT '##INFO=<ID=R_HPA,Number=1,Type=String,Description="RNA expression sum HPA derived">',"\n";
			print OUT '##INFO=<ID=R_GTEx,Number=1,Type=String,Description="RNA expression sum GTEx derived">',"\n";
			print OUT '##INFO=<ID=R_FANTOM,Number=1,Type=String,Description="RNA expression sum FANTOM5 derived">',"\n";
		}
		
		
		print OUT $line;
		#detect the columns that each memeber is in
		chomp $line;
		my @col = split /\t/, $line;
		my( $i )= grep { $col[$_] eq $child } 0..$#col;
		$PR = $i;
		( $i )= grep { $col[$_] eq $mother } 0..$#col;
		$MO = $i;
		( $i )= grep { $col[$_] eq $father } 0..$#col;
		$FR = $i;
	}elsif ($line !~ /^#/){
		$cto++;
		#break up line into information to be used
		chomp $line;
		my @col = split /\t/, $line;
		foreach my $x (0..6){
			print OUT "$col[$x]\t";
		}
		print OUT $col[7];
		my @info = split /;/, $col[7];
		my @format = split /:/, $col[8];
		my( $i )= grep { $format[$_] eq "GT" } 0..$#format;
		my $GT = $i;
		my @cld = split/:/, $col[$PR];
		my @mom = split/:/, $col[$MO];
		my @dad = split/:/, $col[$FR];
		
		if ($freqsw){ #calcualte and annotate frequncies
			foreach my $base(@bases){
				if ($col[7] =~ /\Q$base/){
					my( $i )= grep { $info[$_] =~ /\Q$base/ } 0..$#info;
					my $temp = $info[$i];
					$temp =~ s/\Q$base=//;
					my @temp = split /\//, $temp;
					my @alleles;
					my @counts;
					foreach my $x (0..$#temp){
						my @thing = split /:/, $temp[$x];
						$alleles[$x] = $thing[0];
						$counts[$x] = $thing[1];
					}
					my $total = 0;
					foreach my $number (@counts){
						$total += $number;
					}
					print OUT ";FQ_$base=";
					$freqCt++;
					foreach my $x (0..$#alleles){
						if ($alleles[$x] eq "ref"){
							$alleles[$x] = $col[3];
						}
						print OUT "$alleles[$x]:";
						print OUT ($counts[$x] / $total);
						if ($x < $#alleles){
							print OUT "/";
						}
					}
				}		
			}
		}
		if ($modesw){ #annotate mode of inheirtance
			if ($cld[$GT] !~ /[0-9]/ or $mom[$GT] !~ /[0-9]/ or $dad[$GT] !~ /[0-9]/){
				my @t = err_check($cld[$GT],$mom[$GT],$dad[$GT],$line);
				print OUT ";IM=$t[0]";
				if ($t[0] =~ /^E/){
					$cntE++;
				}elsif ($t[0] =~ /^D/){
					$cntL++;
				}else{
					$cntO++;
				}
			}else{
				my @t = dev_homo_x($cld[$GT],$mom[$GT],$dad[$GT],$male,$col[0]);
				print OUT ";IM=$t[0]_$t[1]";
				if ($t[0] =~ /^D/){
					$cntD++;
				}elsif ($t[0] =~ /^H/){
					$cntR++;
				}elsif ($t[0] =~ /^O/){
					$cntO++;
				}elsif ($t[0] =~ /^X/){
					$cntX++;
				}elsif ($t[0] =~ /^P/){
					$cntP++;
				}elsif ($t[0] =~ /^L/){
					$cntL++;
				}else{
					$cntE++;
				}
			}
		}
		if ($col[7] =~ /gSYM=/ or $col[7] =~ /GL=/){
			my %pgene;
			if ($col[7] =~ /gSYM=([\w\d-]*)/){
				$pgene{$1}=1;
			}if ($col[7] =~ /GL=([\w\d-]*)/){
				$pgene{$1}=1;
			}
			my @pgene = keys %pgene;
			if ($#pgene){
				my @one = split /\t/, $PAT{$pgene[0]};
				my @two = split /\t/, $PAT{$pgene[1]};
				if (defined $one[0] or defined $two[0]){
					print OUT ";PCl=[$pgene[0]/$pgene[1]]$one[0]/$two[0]";
				}if (defined $one[1] or defined $two[1]){
					print OUT ";Pbr=[$pgene[0]/$pgene[1]]$one[1]/$two[1]";
				}if (defined $one[2] or defined $two[2]){
					print OUT ";R_HPA=[$pgene[0]/$pgene[1]]$one[2]/$two[2]";
				}if (defined $one[3] or defined $two[3]){
					print OUT ";R_GTEx=[$pgene[0]/$pgene[1]]$one[3]/$two[3]";
				}if (defined $one[4] or defined $two[4]){
					print OUT ";R_FANTOM=[$pgene[0]/$pgene[1]]$one[4]/$two[4]";
				}
				
			}else{
				my @return = split /\t/, $PAT{$pgene[0]};
				if (defined $return[0]){
					print OUT ";PCl=$return[0]";
				}if (defined $return[1]){
					print OUT ";Pbr=$return[1]";
				}if (defined $return[2]){
					print OUT ";R_HPA=$return[2]";
				}if (defined $return[3]){
					print OUT ";R_GTEx=$return[3]";
				}if (defined $return[4]){
					print OUT ";R_FNATOM=$return[4]";
				}
			}
		}
				
				
		foreach my $x (8..$#col){
			print OUT "\t$col[$x]";
		}
		print OUT "\n";
	}
	if ($screen){
		if ($wirC >= 1){
			print STDERR "Annotating $cto of $entcnt: ", $whirley->(), "\r";
			$wirC = 0;
		}
		else{
			$wirC++;
		}
	}
}
close VAR;
close OUT;


#reoprt
print "\n:::::::::::::::Complete:::::::::::::::\n";
if ($freqsw){
	print "\tAnnotated:\n\t\t$freqCt frequencies\n";
}
if ($modesw){
	print "\tAnnotated:\n\t\t$cntD De novo\n\t\t$cntR Recessive\n\t\t$cntX X-linked\n\t\t$cntP Potentially hemizygous\n\t\t$cntL Potentially Denovo\n\t\t$cntO Other\n\t\t$cntE Error:Null\n";
}


sub err_check{
	if ($_[0] eq './.'){
		return "Er","nul";
	}elsif ($_[1] eq './.'){
		if ($_[0] eq $_[2]){
			return "OT","nul";
		}else{
			my @C = split /\//, $_[0];
			if ($C[0] eq $C[0] and $_[2] !~ /\Q$C[0]/){
				return "Dp", $C[0];
			}
			return "Er","nul";
		}
	}elsif ($_[2] eq './.'){
		if ($_[0] eq $_[1]){
			return "OT","nul";
		}else{
			my @C = split /\//, $_[0];
			if ($C[0] eq $C[0] and $_[1] !~ /\Q$C[0]/){
				return "Dp", $C[0];
			}
			return "Er","nul";
		}
	}else{
		die "\nError bad geneotype: error_check\n@_\n";
	}
}
		

sub dev_homo_x{
	#load in variabls and push all into array
	my $c = $_[0];
	my @C = split /\//, $c;
	if ($_[4] eq "X"){
		#print "X\n";
		if ($_[3] == 1){
			#print "$_[0]\t$_[1]\t$_[2]\n";
			if ($C[0] ne $C[1]){
				return "Er","impossible";
			}else{
				my $M = 0; while($_[1] =~ s/\Q$C[0]//){$M++;}
				my $D = 0; while($_[2] =~ s/\Q$C[0]//){$D++;}
				if ($D == 1){
					return "Er","impossible";
				}elsif ($D == 2){
					return "OT", $C[0];
				}else{
					if ($M == 0){
						return "DX", $C[0];
					}elsif ($M == 1){
						return "XL", $C[0];
					}else{
						return "OT", $C[0];
					}
				}
			}
		}else{
			if ($C[0] eq $C[1] and $C[0] ne "0"){
				my $M = 0; while($_[1] =~ s/\Q$C[0]//){$M++;}
				my $D = 0; while($_[2] =~ s/\Q$C[0]//){$D++;}
				if ($D == 1){
					return "Er","impossible";
				}elsif ($D == 2){
					if ($M == 0){
						return "DN_$C[0],PH", $C[0];
					}else{
						return "OT", $C[0];
					}
				}else{
					return "DN_$C[0],PH", $C[0];
				}
			}else{
				my %return;
				foreach my $C (@C){
					my $M = 0; while($_[1] =~ s/\Q$C//){$M++;}
					my $D = 0; while($_[2] =~ s/\Q$C//){$D++;}
					if ($D == 1){
						return "Er","impossible";
					}elsif ($D == 0 and $M == 0){
						$return{"DE"} = $C;
					}else{
						$return{"OT"} = $C;
					}
				}if (exists $return{"DE"}){
					return "DE", $return{"DE"};
				}else{
					return "OT", $return{"OT"};
				}
			}
		}
	}else{
		if ($C[0] eq $C[1] and $C[0] ne "0"){
			my $M = 0; while($_[1] =~ s/\Q$C[0]//){$M++;}
			my $D = 0; while($_[2] =~ s/\Q$C[0]//){$D++;}
			if ($M == 2 or $D == 2){
				if ($M == 0 or $D == 0){
					return "PH", $C[0];
				}else{
					return "OT", $C[0];
				}
			}elsif ($M == 1 and $D == 1){
				return "HO". $C[0];
			}else{
				return "DN_$C[0],PH", $C[0];
			}
		}else{
			my %return;
			foreach my $C (@C){
				my $M = 0; while($_[1] =~ s/\Q$C//){$M++;}
				my $D = 0; while($_[2] =~ s/\Q$C//){$D++;}
				if ($M == 0 and $D == 0){
					$return{"DE"} = $C;
				}elsif (($M+$D) != 4){
					$return{"OT"} = $C;
				}
			}
			if (exists $return{"DE"}){
				 return "DE", $return{"DE"};
			}else{
				return "OT", $return{"OT"};
			}
		}
	}
}
