#!/usr/bin/perl -w

=head1 LICENSE

Fix_vcfColumn.pl

Copyright (C) 2013 Center for Research Informatics, The University of Chicago

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

############################################################################################
# Riyue Bao 12-04-2012
# This script fixes column title before concatenating vcf files for every chromosome.
############################################################################################


use strict;
use Cwd;
use Cwd 'abs_path';
use IO::Handle;
use Getopt::Long;
use FileHandle;
use File::Basename;


#options
use vars qw($opt_outputDirectory);


##############################
# Arguments from command line
##############################

my $PROG = basename($0);

my $menu = "
#This script removes string \".chr*\" from the column title of every chr*.vcf file before concatenating them together.#
Usage : $PROG -vcf <vcf_file> -o <output_directory>
Input : chr*.vcf
Output : chr*.fixed.vcf

Options :
	[--vcfFile|-vcf] : Input vcf files. More than one vcf file can be specified by multiple -vcf options. 
	[--outputDirectory|-o] : Output directory. Default is the current directory.
	[--verbose|-v] : Print output and error to the screen.
	[--quiet|-q] : Redirect output and error to local log files.
	
Examples : $PROG -vcf \"projPath/multiSampleFiles/projName.chr*.vcf\" -o projPath/multiSampleFile
           $PROG -vcf \"projPath/multiSampleFiles/projName.chr*.mpileup.vcf\" -o projPath/multiSampleFiles
           $PROG -vcf projPath/multiSampleFiles/projName.chr1.GATK.vcf -vcf projPath/multiSampleFiles/projName.chr2.GATK.vcf -o projPath/multiSampleFiles

";

if(@ARGV == 0) { die $menu; }

my $verbose = 1; #print stdout and stderr to the screen or not; turned on by default
my @files;
my $output_dir = "./";

&GetOptions( 
	"verbose|v" => \$verbose,
	"quiet|q" => sub { $verbose = 0 },
	"vcfFile|vcf:s" => \@files,
	"outputDirectory|o:s" => \$output_dir,
	) or die "Unknown option\n";

if (@files == 0) { die "Input vcf files missing. $PROG -vcf <vcf_file> -o <output_directory>\n"; }	

if($opt_outputDirectory) { $output_dir = $opt_outputDirectory; }
if( ! (-d "$output_dir")) { die "Specified output directory does not exist. outputDirectory = $output_dir\n"; }


##############################
# Log files
##############################

my $logFile = $PROG;
$logFile =~ s/.pl$//;
$logFile = "$logFile.log";
if($verbose)
{
	#print stdout and stderr to the screen
}
else
{
	open(STDOUT, ">", "$logFile.stdout") or die $!;
	open(STDERR, ">", "$logFile.stderr") or die $!;
}

##############################
# Main body
##############################

print "\n[COMMANDS]\n$PROG --vcfFile ";
foreach my $file (@files) { print "$file "; }
print "--outputDirectory $output_dir ";
if($verbose) { print "--verbose " }
else { print "--quiet "; }
print "\n";
print "\n[PROGRESS]\n";

$output_dir =~ s/\/$//;
if(@files == 1 && $files[0] =~ m/^(\S*)\*(\S*)$/) { @files = <$1*$2>; }

Update_Progress("Program run starts");

my $file_total = 0;
my $outputfile_total = 0;
foreach my $file (@files) {
	if($file !~ m/fixed/)
	{
		$file_total++;
		my $output_file = $file;
		if($output_file =~ m/^\S+\//) { $output_file =~ s/^\S+\///g; }
		$output_file =~ s/vcf$/fixed.vcf/;
		$output_file = "$output_dir/$output_file";
		
		Update_Progress("Processing \"$file\"; Output written into \"$output_file\"");
		$file = Open_File($file, $output_file);
	}
}

Update_Progress("$file_total input vcf files found; $outputfile_total output fixed.vcf files written");
if($file_total != $outputfile_total) { print STDERR "Total number of input files is not equal to that of output files. Failed to write some output files?\n"; }
Update_Progress("Program run finished\n");

##############################
# Functional Modules
##############################

sub Update_Progress
{
	my $progress = shift;
	
	print "[",Local_Time(), "]$progress\n";
	
	return $progress;
}

sub Open_File
{
	my ($file, $output_file) = @_;
	
	open(IN, "<", $file) or die $!;
	open(OUT, ">", $output_file) or die $!;
	
	while(my $line = <IN>) {
		chomp($line);
		
		if($line =~ m/^#CHROM/)
		{
			my @line = split(/\t/, $line);
			for(my $i = 9; $i < @line; $i++) {
				if($line[$i] =~ m/^(\S+)\.chr/)
				{
					$line[$i] = $1;
				}
				if($line[$i] =~ m/^\S+\/\S+/)
				{
					my @temp = split(/\//, $line[$i]);
					$line[$i] = $temp[@temp-2];
				}
			
			}
			$line = join("\t", @line);
			
		}
		elsif($line =~ m/^[1-9XYM]+/)
		{
			my @line = split(/\t/, $line);
			$line[0] = "chr".$line[0];
			$line = join("\t", @line);
		}
		
		print OUT "$line\n";
	
	}
	
	close(IN);
	close(OUT);
	
	$outputfile_total++;
	
	return $file;
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

