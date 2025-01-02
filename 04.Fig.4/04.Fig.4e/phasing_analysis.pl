#!/usr/bin/perl
#use strict;
use Benchmark;
# modified by FL 20190317
#----------------------------------------------------
my $abundance = $ARGV[0];
my $phase = $ARGV[1];

my $cycles = 10; #define the cycles for caculation
my $minimal_n = 3;
my $sum = 0;

#----------------------------------------------------
open(ABUN, $abundance) or die("Could not open  file.");
while(my $line = <ABUN>)  {  
    chomp($line);
    my ($chr,$position,$abundance) = split(/\s/,$line);
    $sum += $abundance;
    my @list = ($chr, $position, $abundance);
    push(@all_position_abun,[@list]);
}
undef %hash_formated_position;
@all_position_abun_sorted = sort {$a->[0] <=> $b->[0] or $a->[1] <=> $b->[1]} @all_position_abun;
#print "\@all_position_abun sorted\n";
#----------------------------------------------------

    open (PHASED_PST, ">./phase\_score.txt");
    open (ALL_PST, ">./phase\_All\_score.txt");
    print PHASED_PST "Chr\tStart\tAbundance\tWindow_N_loci\tPhased_N_position\tPhased_N_loci\tTotal_count\tPhased_count\tMax_phased_count\tPhased_ratio\tMax_phased_ratio\tPhasing_Score\n";
    print ALL_PST "Chr\tStart\tAbundance\tWindow_N_loci\tPhased_N_position\tPhased_N_loci\tTotal_count\tPhased_count\tMax_phased_count\tPhased_ratio\tMax_phased_ratio\tPhasing_Score\n";
    my $chr_current = 1; #initialize the current chromosome number
    my @window = (); #initialize the calculation window
    foreach my $key (@all_position_abun_sorted){
		my @array = @{$key};
        	$array[0] =~ s/$array[0]/1/i;
		my $chr = $array[0];
		my $pst = $array[1];
		if (@window == 0){
		push(@window,[@array]);
		next;
		}
		$first_pst = $window[0][1];
		$window_end = $first_pst + $cycles*$phase - 1;
		if ($chr == $chr_current){	
			if ($pst <= $window_end) { #read smallRNA into the window
				push(@window,[@array]);
			}
			else {
				#select the phased sRNA reads into @phased matrix
				for (my $i = 0; $i < @window; $i++) {
					if (($window[$i][1]-$first_pst)%$phase == 0) { 
					push (@phased, [$window[$i][0], $window[$i][1],$window[$i][2]]);
					$hash_3{$window[$i][1]} = $array[2];
					}
				 	elsif (($window[$i][1] + 1 -$first_pst)%$phase == 0) { 
				        push (@phased, [$window[$i][0], $window[$i][1],$window[$i][2]]);
			        	$hash_3{$window[$i][1] + 1} = $array[2];
					}
					elsif (($window[$i][1] - 1 -$first_pst)%$phase == 0) { 
			        	push (@phased, [$window[$i][0], $window[$i][1],$window[$i][2]]);
					$hash_3{$window[$i][1] - 1} = $array[2];
					}
				}
				#sum the abundance for phased positions
				for (my $i = 0; $i< @phased ; $i++) {
					if ($phased[$i][-1] >= 1){
						$k +=1;
						$sum_ki += $phased[$i][-1];
       					}
					if ($kp_max <= $phased[$i][-1]) {
						$kp_max = $phased[$i][-1];
					}
				}
				#sum the abundance of total positions
				for (my $j = 0; $j< @window ; $j++) {
					$sum_total += $window[$j][-1];
				}
				$n = @window;
				$k_p = keys %hash_3;    #number of phased loci (with tolerance to -1 and +1)
				if ($k_p >= $minimal_n ) { 
					$P_value = log((1+($sum_ki))**($k-2));
					$P_value = sprintf "%.2f", $P_value;
				}
				else {
					$P_value = 0;
				}
				$output = join ("\t", @{$window[0]});
				$phased_ratio = $sum_ki/$sum;
				if ($sum_ki > 0 ){$max_ratio = $kp_max/$sum_ki;}else{ $max_ratio = NA;}
				print PHASED_PST "$output\t$n\t$k_p\t$k\t$sum_total\t$sum_ki\t$kp_max\t$phased_ratio\t$max_ratio\t$P_value\n" if ($P_value >0);
				print ALL_PST "$output\t$n\t$k_p\t$k\t$sum_total\t$sum_ki\t$kp_max\t$phased_ratio\t$max_ratio\t$P_value\n";
				shift @window;
				@phased = ();
				$sum_ki = 0;
				$sum_total = 0;
				$P_value= 0;
				$n = 0;
				$k = 0;
				$k_p = 0;
				%hash_3 = ();
				$kp_max = 0;
				redo;
			}
		}
		else {
		$chr_current++;
		@window = ();
		push(@window,[@array]);
		} 
	
	}
	close IN3;
	close PHASED_PST;
	close ALL_PST;
	close ALL_PST;
