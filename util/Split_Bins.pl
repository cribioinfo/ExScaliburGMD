#!/usr/bin/perl -w

=head1 LICENSE

Split_Bins.pl

Copyright (C) 2013-2015 Center for Research Informatics, The University of Chicago

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

############################################################################################
# Riyue Bao 01/26/2015
# This script takes a sorted bed file, concats the regions and splits by even-sized bins. 
############################################################################################


use strict;
use Cwd;
use Cwd 'abs_path';
use IO::Handle;
use Getopt::Long;
use FileHandle;
use File::Basename;

##############################
# Arguments from command line
##############################

my $PROG = basename($0);

my $menu = qq~
## This script takes a bed file, concats the regions and splits by even-sized bins. 
## The bed file input must be sorted! (e.g. use bedtools sort) 
Usage: 
$PROG -f|--file <sorted_bedfile> -o|--output <outputfile> -n|--bin <bin> 

Example: 
$PROG -f LCAexome.bwamem.target.bed -o LCAexome.bwamem.target.bins.txt -b 25000

~;

if(@ARGV == 0) { print "\n$menu\n"; exit; }

my $file = "";
my $output_file = "bins.txt";
my $bin = 1000;

&GetOptions( 
	"file|f:s" => \$file,
	"output|o:s" => \$output_file,
	"bin|b:i" => \$bin
) or die $!;

if ($file eq "") { print STDERR "Input file missing! \n$menu\n"; exit(1); }	

my $command = "$PROG --file $file --output $output_file --bin $bin";


##############################
# Main 
##############################

print"\n[COMMANDS]\n$command\n\n[PROGRESS]\n";

Update_Progress("Opening input bed file $file");
$output_file = Open_File($file, $output_file);

Update_Progress("Output file written into $output_file");

Update_Progress("Program run finished\n");

##############################
# Functions
##############################

sub Open_File
{
	my ($file, $output_file) = @_;

	open(IN, "<", $file) or die $!;
	open(OUT, ">", $output_file) or die $!;

	print OUT "#CHROM\tSTART\tEND\tBIN\n";

	my $i = 0;
	my ($chr_prev, $start_prev, $end_prev);
	my $chr_regions;
	my $chr_counts;
	my $chr_bins;
	my $bin_count = 0;
	my $bin_length = 0; ## keep track of how long the region is coverted from the first region

	while(my $line = <IN>) {
		chomp($line);
		my @line = split(/\t/, $line);

		if($line =~ m/^chr/) {
			$i++;

			my ($chr, $start, $end) = @line;
			$chr_regions->{$chr} += $end - $start;
			$chr_counts->{$chr}++;

			if($i == 1) {
				$bin_length = 0;
				$bin_count = 0;
				$chr_bins->{$chr}++;
				($bin_length, $bin_count, $chr_bins) = Print_Bin($bin_length, $bin_count, $chr, $start, $end, $line, $chr_bins);

			}
			if($i > 1){
				if($chr eq $chr_prev) {
					
					($bin_length, $bin_count, $chr_bins) = Print_Bin($bin_length, $bin_count, $chr, $start, $end, $line, $chr_bins);
					
				}
				else {
					## start a new chrom
					$bin_length = 0;
					$bin_count = 0;
					$chr_bins->{$chr}++;
					($bin_length, $bin_count, $chr_bins) = Print_Bin($bin_length, $bin_count, $chr, $start, $end, $line, $chr_bins);
				}
			}

			($chr_prev, $start_prev, $end_prev) = ($chr, $start, $end);

		}

	}

	close(IN);
	close(OUT);

	print "lineTotal = $i\n";
	print"#Chr\tRegionCount\tRegionLength\tBinCount(bin=$bin"."bp)\n";
	foreach my $chr(sort keys %{$chr_counts}) {
		print "$chr\t".$chr_counts->{$chr};
		if(exists $chr_regions->{$chr}) {
			print"\t".$chr_regions->{$chr};
		}
		if(exists $chr_bins->{$chr}) {
			print"\t".$chr_bins->{$chr};
		}
		print "\n";
	}

	return $output_file;
}

sub Print_Bin
{
	my ($bin_length, $bin_count, $chr, $start, $end, $line, $chr_bins) = @_;

	$bin_length += $end - $start;
	# print "$bin_count : $bin_length\n";
	if($bin_length < $bin * ($bin_count+1)) {
		print OUT "$line\t$bin_count\n";
	}
	else {
		my $count1 = int(($bin_length - ($end - $start)) / $bin);
		my $count2 = int($bin_length / $bin);
		if((($bin_length - ($end - $start)) % $bin) > 0 ) { $count1 += 1; }
		# print "$count1, $count2\n";

		my $start1 = $start;
		my $end1 = $start1 + ($bin * $count1 - ($bin_length - ($end - $start)));
		
		## print the head 
		if($start1 != $end1) {
			print OUT "$chr\t$start1\t$end1\t$bin_count\n";
		}

		## print the middle 
		for(my $j=0; $j<$count2-$count1;$j++) {
			$bin_count++;
			$start1 = $end1;
			$end1 = $start1 + $bin;
			$chr_bins->{$chr}++;

			print OUT "$chr\t$start1\t$end1\t$bin_count\n";
		}

		## print the tail 
		if($end1 < $end) {
			$bin_count++;
			$chr_bins->{$chr}++;
			print OUT "$chr\t$end1\t$end\t$bin_count\n"; 

		}
		elsif($end1 > $end) 
		{
			print STDERR "$line : end1 = $end1, end = $end\n";
		}

	}

	return $bin_length, $bin_count, $chr_bins;
}


sub Update_Progress
{
	my $progress = shift;
	
	print "[",Local_Time(), "]$progress\n";
	
	return $progress;
}


sub Test_Hash
{
	my $hash = shift;
	
	foreach my $key (sort keys %{$hash}) {
		print "$key=>\n";
		
		if($hash->{$key} =~ m/^HASH/)
		{
			Test_Hash($hash->{$key});
			print "\n";
		}
		elsif($hash->{$key} =~ m/^ARRAY/)
		{
			print join(" ", @{$hash->{$key}}), "\n";
		}
		else
		{
			print $hash->{$key}, "\n";
		}
	
	}

	return $hash;
}

#Print local time
sub Local_Time
{
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my  @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	#my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek], $months[$month] $dayOfMonth, $year";
	
	my $digits = 2;
	$month = Add_Prefix($month, $digits);
	$dayOfMonth = Add_Prefix($dayOfMonth, $digits);
	$hour = Add_Prefix($hour, $digits);
	$minute = Add_Prefix($minute, $digits);
	$second = Add_Prefix($second, $digits);	
	my $theTime = "$month-$dayOfMonth-$year\_$hour:$minute:$second";
	
	#print "Local time: ", $theTime, "\n";
	return $theTime;
} 

sub Add_Prefix
{
	my ($number, $digits) = @_;
	
	if($number =~ m/^\d{$digits}$/)
	{	
		return $number;
	}
	else
	{
		$number = "0".$number;
		Add_Prefix($number, $digits);
	}

}

