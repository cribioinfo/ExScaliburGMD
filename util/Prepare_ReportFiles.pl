#!/usr/bin/perl -w

=head1 LICENSE

Prepare_ReportFiles.pl

Copyright (C) 2013-2015 Center for Research Informatics, The University of Chicago

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

############################################################################################
# Riyue Bao 02/26/2015
# This script prepares and processes fastqc/aln/annovar files as input for pipeline report utility.
############################################################################################


use strict;
use Cwd;
use Cwd 'abs_path';
use IO::Handle;
use Getopt::Long;
use FileHandle;
use File::Basename;
use YAML::Tiny;
use Data::Dumper;

## options
use vars qw($opt_pass);

##############################
# Arguments from command line
##############################

## print menu
my $PROG = basename($0);
my $menu = Print_Menu($PROG);
if(@ARGV == 0) { print "\n$menu\n"; exit; }

## get opt 
my $in = "";
my $out = "";
my $outputdir = "./";
my $pass_flag = 0;

&GetOptions( 
	"in|i:s" => \$in,
	"out|o:s" => \$out,
	"dir|d:s" => \$outputdir,
	"pass|p"
) or die $!;

## check
if ($in eq "" || $out eq "") { print STDERR "Input file missing! \n$menu\n"; exit(1); }	
if($opt_pass) { $pass_flag = 1; }

## command
my $command = "$PROG --in $in --out $out --dir $outputdir";
if($pass_flag) { $command .= " --pass"; }

$outputdir =~ s/\/$//;
if ( ! -d $outputdir) { mkdir $outputdir; }

##############################
# Main 
##############################

print "\n[COMMANDS]\n$command\n\n[PROGRESS]\n";

## initialize 
my $in_list;
my $out_list;
my $variant_list; ## combine various variant sets for printing
my $report_list; ## read in report output file info for printing
my $sample_list; 
my @aligners = qw(bwaaln bwamem novoalign);
my @callers = qw(gatkhc gatkug freebayes mpileup ivc platypus vcmm);

# print Dumper($in_list);
# print Dumper($out_list);
# print $in_list->{"project"}."\n";
# YAML::Tiny::DumpFile("test.out",$in_list);

## run start
Update_Progress("Read input file info from $in");
$in_list = YAML::Tiny::LoadFile($in);

Update_Progress("Read output file info from $out");
$out_list = YAML::Tiny::LoadFile($out);

## read in report output file info
Update_Progress("Prepare sample2file list");
$report_list = Parse_ReportOut($out_list, $report_list);
# Test_Hash($report_list);

## fastqc
Update_Progress("Parse fastqc files for QC stats");
$in_list = Parse_FastQC($in_list, $report_list);

## metrics
Update_Progress("Parse alignment metrics and coverage");
$in_list = Parse_Metrics($in_list, $report_list);

## annovar
Update_Progress("Parse variant files for multiple aligners/callers");
$in_list = Parse_Annovar($in_list, $report_list, $variant_list, $sample_list);

## run end
Update_Progress("Output files written into directory $outputdir");
Update_Progress("Program run finished\n");

##############################
# Functions
##############################

sub Parse_Metrics
{
	my ($in_list, $report_list) = @_;

	my @array1 = qw( TOTAL_READS PF_READS_ALIGNED PCT_PF_READS_ALIGNED PF_HQ_ALIGNED_READS PF_MISMATCH_RATE PF_HQ_ERROR_RATE PF_INDEL_RATE MEAN_READ_LENGTH PCT_READS_ALIGNED_IN_PAIRS STRAND_BALANCE PCT_CHIMERAS );
	my @array2 = qw(MEDIAN_INSERT_SIZE MEDIAN_ABSOLUTE_DEVIATION MIN_INSERT_SIZE MAX_INSERT_SIZE MEAN_INSERT_SIZE STANDARD_DEVIATION READ_PAIRS PAIR_ORIENTATION WIDTH_OF_10_PERCENT WIDTH_OF_20_PERCENT WIDTH_OF_30_PERCENT WIDTH_OF_40_PERCENT WIDTH_OF_50_PERCENT WIDTH_OF_60_PERCENT WIDTH_OF_70_PERCENT WIDTH_OF_80_PERCENT WIDTH_OF_90_PERCENT WIDTH_OF_99_PERCENT);
	my @array3 = qw();
	my @array4 = qw();
	my @array5 = qw();

	my $metrics_list = { 
		"alignment_summary_metrics" => \@array1,
		"insert_size_metrics" => \@array2,
		"quality_by_cycle_metrics" => \@array3,
		"quality_distribution_metrics" => \@array4,
		"total_coverage" => \@array5
		};
	# print "Metrics = ".join(" ", keys %{$metrics_list})."\n";

	if(exists $in_list->{"data"}) {
		my @data = @{$in_list->{"data"}};
		
		foreach my $data (@data) { 
			# print $data[2]->{"variants"}->{"inputs"}->{"variant_table"}."\n";
			my $sample = "";
			if(exists $data->{"sample"}) { $sample = $data->{"sample"}; }
			print "-----------\n";
			print "sample = $sample\n";

			if($sample =~ m/\S+/) {
				if(exists $data->{"alignments"}->{"inputs"}) {
					my $path = "";
					my @files;
					my @keys  = qw(alignment_summary_metrics insert_size_metrics quality_by_cycle_metrics quality_distribution_metrics total_coverage);
					if(exists $data->{"alignments"}->{"inputs"}->{"path"}) { 
						$path = $data->{"alignments"}->{"inputs"}->{"path"}; 
					}
					if($path !~ m/\S+/) {
						print STDERR "Parse_Metrics: in_list:data:alignments:inputs:path: key does not exist!\n";
					}
					else 
					{
						foreach my $key (@keys) {
							if(exists $data->{"alignments"}->{"inputs"}->{$key}) { 
								@files = split(/,/, $data->{"alignments"}->{"inputs"}->{$key});
								for(my $i=0; $i<@files; $i++) { $files[$i] = "$path/$files[$i]"; }
								if(exists $report_list->{$sample}->{"alignments"}->{$sample}->{$key}) {
									my $printfile = "$outputdir/".$report_list->{$sample}->{"alignments"}->{$sample}->{$key};
									$printfile = Print_Metrics($printfile, $metrics_list, $key, $sample, \@files);
								}
								else {
									print STDERR "Parse_Metrics:report_list:$sample:alignments:$sample:$key: key does not exist!\n";
								}
							}
							else {
								print STDERR "Parse_Metrics:in_list:alignments:inputs:$key: key does not exist!\n";
							}
						}
					}
					
				}
				else {
					print STDERR "Parse_Metrics: in_list:data:alignments:inputs: key does not exist!\n";
				}
			}
			else { print STDERR "Parse_Metrics: in_list:data:alignments:sample: value is empty!\n";}

		}
	}
	else {
		print STDERR "Parse_Metrics:in_list:data: key does not exist!\n";
	}


	return $in_list;
}

sub Print_Metrics
{
	my ($printfile, $metrics_list, $metrics, $sample, $files_ref) = @_;

	my @files = @{$files_ref};
	my $type = "NORMAL";
	my $printdir = "";
	my $i=0;

	if($printfile =~ m/(^\S+\/)/g) { $printdir = $1; }
	if( ! -d $printdir) { `mkdir -p $printdir`; }

	open(OUT, ">", $printfile) or die $!;
	foreach my $aligner (@aligners) {
		foreach my $file(@files) {
			if($file =~ m/\.$aligner\.\S*$metrics/i || $file =~ m/\.$aligner\.\S*covhist/i) {
				$i++;
				my @categories = @{$metrics_list->{$metrics}};
				print "-----------\n";
				print "metrics = $metrics\naligner = $aligner\n";
				print "file = $file\n";
				$file = Open_File_Metrics($file, $metrics, $aligner, $sample, $type, $i, \@categories);

			}
		}
	}
	close(OUT);
	print "Output written into $printfile\n";

	return $printfile;
}

sub Open_File_Metrics
{
	my ($file, $metrics, $aligner, $sample, $type, $i, $categories_ref) = @_;
	
	open(IN, "<", $file) or die $!;

	my @header;
	my @categories = @{$categories_ref};
	# print "Categories include:\n".join("\n", @categories)."\n";
	## notice that array 3 and 4 are empty so you won't have any categories printed out.

	if($metrics eq "alignment_summary_metrics")
	{
		@header = qw(Sample Type Aln Metric Category Measure Value);
		if($i==1) { print OUT join("\t", @header)."\n"; }
		
		my $field_list;
		my $key_list = {
			"FIRST_OF_PAIR" => "R1",
			"SECOND_OF_PAIR" => "R2",
			"PAIR" => "Paired",
			"UNPAIRED" => "Unpaired"
		
		};
		
		while(my $line = <IN>) {
			chomp($line);
			my @line = split(/\t/, $line);
			
			if($line =~ m/^CATEGORY/)
			{
				## only include categories predefined in @categories
				for(my $i=0; $i<@line; $i++) {
					my $field = $line[$i];
					if(grep(/^$field$/, @categories)) ## if those category is from @categories
					{
						$field_list->{$i} = $field;
						# print "$i => $field\n";
					}
				}
			}
			
			foreach my $key (sort keys %{$key_list}) {
				my $value = $key_list->{$key};
				if($line =~ m/^$key/)
				{
					foreach my $i (sort{$a<=>$b} keys %{$field_list}) {
						print OUT "$sample\t$type\t$aligner\t$metrics\t$value\t".$field_list->{$i}."\t$line[$i]\n";
					}
				}
			}
	
		}
	}
	elsif($metrics eq "insert_size_metrics")
	{
		@header = qw(Sample Type Aln Metric Category Measure Value);
		if($i==1) { print OUT join("\t", @header)."\n"; }
		
		my $field_list;
		my $table_line = "";
		my $hist_line = "";
		my $i = 0;
		
		while(my $line = <IN>) {
			chomp($line);
			my @line = split(/\t/, $line);
			$i++;
			
			if($line =~ m/^MEDIAN_INSERT_SIZE/)
			{
				$table_line = $i;
				## only include categories predefined in @categories
				for(my $i=0; $i<@line; $i++) {
					my $field = $line[$i];
					if(grep(/^$field$/, @categories)) ## if those category is from @categories
					{
						$field_list->{$i} = $field;
						# print "$i => $field\n";
					}
				}
			}
			
			if($table_line ne "" && $i == $table_line + 1)
			{
				foreach my $i (sort{$a<=>$b} keys %{$field_list}) {
					print OUT "$sample\t$type\t$aligner\t$metrics\tTABLE\t".$field_list->{$i}."\t$line[$i]\n";
				}
			}

			if($line =~ m/^insert_size/)
			{
				$hist_line = $i;
			}
			
			if($hist_line ne "" && $i >= $hist_line + 1 && $line =~ m/^\d+/)
			{
				print OUT "$sample\t$type\t$aligner\t$metrics\tHIST\t$line\n";
			}			
		}
	}
	elsif($metrics eq "quality_by_cycle_metrics")
	{
		@header = qw(Sample Type Aln Metric Cycle Quality);
		if($i==1) { print OUT join("\t", @header)."\n"; }
		
		my $field_list;
		my $quality_line = "";
		my $i = 0;
		
		while(my $line = <IN>) {
			chomp($line);
			my @line = split(/\t/, $line);
			$i++;


			if($line =~ m/^CYCLE/)
			{
				$quality_line = $i;
			}
			
			if($quality_line ne "" && $i >= $quality_line + 1 && $line =~ m/^\d+/)
			{
				print OUT "$sample\t$type\t$aligner\t$metrics\t$line\n";
			}			
		}
	}
	elsif($metrics eq "quality_distribution_metrics")
	{
		@header = qw(Sample Type Aln Metric Quality Count);
		if($i==1) { print OUT join("\t", @header)."\n"; }
		
		my $field_list;
		my $quality_line = "";
		my $i = 0;
		
		while(my $line = <IN>) {
			chomp($line);
			my @line = split(/\t/, $line);
			$i++;


			if($line =~ m/^QUALITY/)
			{
				$quality_line = $i;
			}
			
			if($quality_line ne "" && $i >= $quality_line + 1 && $line =~ m/^\d+/)
			{
				print OUT "$sample\t$type\t$aligner\t$metrics\t$line\n";
			}			
		}
	}	
	elsif($metrics eq "total_coverage")
	{
		@header = qw(Sample Type Aln Coverage);
		if($i==1) { print OUT join("\t", @header)."\n"; }
		
		my $field_list;
		my $total_length = 0;
		my $total_cov = 0;
		my $i = 0;
		
		while(my $line = <IN>) {
			chomp($line);
			my @line = split(/\t/, $line);
			$i++;

			if($line =~ m/^all/)
			{
				$total_length += $line[2];
				$total_cov += $line[1] * $line[2];
			}		
		}

		print OUT "$sample\t$type\t$aligner\t".sprintf("%.2f",($total_cov/$total_length))."\n";
	}		
	
	close(IN);
	return $file;
}

sub Parse_FastQC
{
	my ($in_list, $report_list) = @_;

	if(exists $in_list->{"data"}) {
		my @data = @{$in_list->{"data"}};
		
		foreach my $data (@data) { 
			# print $data[2]->{"variants"}->{"inputs"}->{"variant_table"}."\n";
			my $sample = "";
			if(exists $data->{"sample"}) { $sample = $data->{"sample"}; }
			print "-----------\n";
			print "sample = $sample\n";

			if($sample =~ m/\S+/) {
				if(exists $data->{"fastqc"}->{"inputs"}) {
					my @inputs = @{$data->{"fastqc"}->{"inputs"}};
					foreach my $input (@inputs) {
						my $rg = "";
						my $path = "";
						my @keys = ("leftseq", "rightseq");
						if(exists $input->{"readgroup"}) { $rg = $input->{"readgroup"}; }
						if(exists $input->{"path"}) { $path = $input->{"path"}; }
						print "-----------\n";
						print "readgroup = $rg\n";

						if($rg !~ m/\S+/) { print STDERR "Parse_FastQC:in_list:data:fastqc:inputs:readgroup: value is empty!" }
						if($path !~ m/\S+/) { print STDERR "Parse_FastQC:in_list:data:fastqc:inputs:path: value is empty!" }
						if($rg =~ m/\S+/ && $path =~ m/\S+/) {
							foreach my $key (@keys) {
								my $file = "";
								if(exists $input->{$key}) { 
									$file = "$path/".$input->{$key}."/fastqc_data.txt"; 
									if(exists $report_list->{$sample}->{"fastqc"}->{$rg}->{$key}) {
										my $printdir = "$outputdir/".$report_list->{$sample}->{"fastqc"}->{$rg}->{$key};
										$printdir = Print_FastQC($printdir, $file);
									}
									else {
										print STDERR "Parse_FastQC:report_list:$sample:fastqc:$rg:$key: key does not exist!\n";
									}
								}
							}
						}
						
					}
					
				}
				else {
					print STDERR "Parse_FastQC: in_list:data:fastqc:inputs: key does not exist!\n";
				}
			}
			else { print STDERR "Parse_FastQC: in_list:data:sample: value is empty!\n";}

		}
	}
	else {
		print STDERR "Parse_FastQC:in_list:data: key does not exist!\n";
	}

	return $in_list;
}

sub Print_FastQC
{
	my ($printdir, $file) = @_;

	if ( ! -d $printdir) { `mkdir -p $printdir` ; }
	print "file = $file\n";

	my @array1 = qw( Base Mean Median LQ UQ P10 P90 );
	my @array2 = qw( Base Count );
	# my @array3 = qw( Base PGC ); ## with fastqc/0.10.1 (older version)
	my @array3 = qw( PGC Count ); ## with fastqc/0.11.2 (the GC title changed!)
	my @array4 = qw();

	my $metrics_list = { 
		"fastqc_qbd" => \@array1,
		"fastqc_pbnc" => \@array2,
		"fastqc_gcbd" => \@array3,
		"fastqc_data" => \@array4,
		};
	print "FastQC stats = ".join(" ", keys %{$metrics_list})."\n";


	## print fastqc report output files
	open(IN, "<", $file) or die $!;

	my @header;
	#my @categories = @$categories_ref;
	#print "Categories include:\n".join("\n", @categories)."\n";

	my $fh_list;
	print "Output written into: \n";
	foreach my $metrics (sort keys %{$metrics_list}) {
		my $output_file = "$printdir/$metrics.txt";
		print "$output_file\n";
		open(my $fh, ">", $output_file) or die $!;
		$fh_list->{$metrics} = $fh;
	}
	#Test_Hash($fh_list);
	
	my $fh;
	my $metrics = "";
	while(my $line = <IN>) {
		chomp($line);
		
		if($line =~ m/^##FastQC/ && $metrics eq "")
		{
			$metrics = "fastqc_data";
			$fh = $fh_list->{$metrics};
		}	
		elsif($line =~ m/^\>\>END_MODULE/ && $metrics eq "fastqc_data")
		{
			$metrics = "";
		}
		elsif($line =~ m/^#Base\s+Mean\s+Median\s+Lower\s+Quartile\s+Upper\s+Quartile/ && $metrics eq "")
		{
			$metrics = "fastqc_qbd";
			$fh = $fh_list->{$metrics};
			my @header = @{$metrics_list->{$metrics}};
			print $fh join("\t", @header)."\n";
		}
		elsif($line =~ m/^\d+/ && $metrics eq "fastqc_qbd")
		{
			print $fh "$line\n";
		}
		elsif($line =~ m/^\>\>END_MODULE/ && $metrics eq "fastqc_qbd")
		{
			$metrics = "";
		}
		# elsif($line =~ m/^#Base\s+\%GC/ && $metrics eq "")
		elsif($line =~ m/^#GC\s+Content\s+Count/ && $metrics eq "")
		{
			print "111111\n";			
			$metrics = "fastqc_gcbd";
			$fh = $fh_list->{$metrics};
			my @header = @{$metrics_list->{$metrics}};
			print $fh join("\t", @header)."\n";
		}
		elsif($line =~ m/^\d+/ && $metrics eq "fastqc_gcbd")
		{
			print $fh "$line\n";
		}
		elsif($line =~ m/^\>\>END_MODULE/ && $metrics eq "fastqc_gcbd")
		{
			$metrics = "";
		}
		elsif($line =~ m/^#Base\s+N\-Count/ && $metrics eq "")
		{
			$metrics = "fastqc_pbnc";
			$fh = $fh_list->{$metrics};
			my @header = @{$metrics_list->{$metrics}};
			print $fh join("\t", @header)."\n";
		}
		elsif($line =~ m/^\d+/ && $metrics eq "fastqc_pbnc")
		{
			print $fh "$line\n";
		}
		elsif($line =~ m/^\>\>END_MODULE/ && $metrics eq "fastqc_pbnc")
		{
			$metrics = "";
		}
		
		if($line =~ m/^\S+/ && $metrics eq "fastqc_data")
		{
			print $fh "$line\n";
		}	
	}
	

	close(IN);
	
	foreach my $metrics (sort keys %{$metrics_list}) {
		my $output_file = "$printdir/$metrics.txt";
		my $fh = $fh_list->{$metrics};
		close($fh);
	}	

	return $file;
}

sub Parse_Annovar
{
	my ($in_list, $report_list, $variant_list, $sample_list) = @_;

	if(exists $in_list->{"multi"}->{"variants"}->{"inputs"}) {
		my @inputs = @{$in_list->{"multi"}->{"variants"}->{"inputs"}};
		my @files;
		my $algorithm_list; ## record which aligner/caller was included in analysis

		foreach my $input (@inputs) {
			my $path = "";
			my $file = "";
			if(exists $input->{"path"}) { $path = $input->{"path"}; }
			if(exists $input->{"variant_table"}) { $file = $input->{"variant_table"}; }

			if($path =~ m/\S+/ && $file =~ m/\S+/) { push(@files, "$path/$file"); }
			else { print STDERR "Parse_Annovar:in_list:multi:variants:inputs:variant_table: value is empty!\n"; }

		}

		print scalar(@files)." annotated variant sets found!\n";
		# print join("\n", @files)."\n";
		foreach my $aligner (@aligners) {
			foreach my $caller (@callers) {
				foreach my $file (@files) {
					if($file =~ m/\.$aligner\.$caller\./i) {
						$algorithm_list->{$aligner}->{$caller} = $file;
						print "-----------\n";
						print "aligner = $aligner\ncaller = $caller\nfile = $file\n";
						($variant_list, $sample_list) = Open_File_Annovar($variant_list, $sample_list, $file, $aligner, $caller);
					}
				}
			}
		}
		# Test_Hash($variant_list);
		# Test_Hash($sample_list);

		print "Print variant per sample ...\n";
		$variant_list = Print_Annovar($variant_list, $sample_list, $algorithm_list,  $report_list);

	}
	else {
		print STDERR "Parse_Annovar:in_list:multi:variants:inputs: key does not exist!\n";
	}

	return $in_list;
}

## refGene, snp138, esp6500siv2_all, 1000g2014sep, cosmic70, clinvar_20140929
sub Print_Annovar
{
	my ($variant_list, $sample_list, $algorithm_list, $report_list) = @_;

	my $field_list = Annovar_Field();
	# my @fields = qw(Chr Start End Ref Alt Func Gene ExonicFunc AAChange dbSNP ESP6500 Genome1000 ExAC65000 COSMIC CADD);
	my @fields = qw(Chr Start End Ref Alt Func Gene ExonicFunc AAChange dbSNP ESP6500 Genome1000 ExAC65000 COSMIC CADD ClinVar);
	my @header = ("Sample", @fields, qw(TotalScore Aligner Caller)); 
	foreach my $aligner (sort keys %{$algorithm_list}) {
		foreach my $caller (sort keys %{$algorithm_list->{$aligner}}) {
			push(@header, "$aligner.$caller");
		}
	}

	foreach my $sample (sort keys %{$sample_list}) {
		my $key = "variant_table";
		if(exists $report_list->{$sample}->{"variants"}->{$sample}->{$key}) {
			my $printfile = "$outputdir/".$report_list->{$sample}->{"variants"}->{$sample}->{$key};
			my $printdir;
			if($printfile =~ m/^(\S+\/)/g) { $printdir = $1; }
			if( ! -d $printdir) { `mkdir -p $printdir`; }
			
			open(OUT, ">", $printfile) or die $!;
			print OUT join("\t", @header)."\n";
			foreach my $chr (sort keys %{$sample_list->{$sample}}) {
				foreach my $start (sort{$a<=>$b} keys %{$sample_list->{$sample}->{$chr}}) {
					foreach my $end (sort{$a<=>$b} keys %{$sample_list->{$sample}->{$chr}->{$start}}) {
						foreach my $ref (sort keys %{$sample_list->{$sample}->{$chr}->{$start}->{$end}}) {
							foreach my $alt (sort keys %{$sample_list->{$sample}->{$chr}->{$start}->{$end}->{$ref}}) {
								my @print = ($sample);
								my $total_score = 0;
								my $aligner_list;
								my $caller_list;

								## print variant annotation fields
								foreach my $field (@fields) {
									if(exists $field_list->{$field}) {
										my $key = $field_list->{$field};
										if(exists $variant_list->{$chr}->{$start}->{$end}->{$ref}->{$alt}->{$key}) {
											my $value = $variant_list->{$chr}->{$start}->{$end}->{$ref}->{$alt}->{$key};
											if($value !~ m/\S+/) { $value = "NA"; }
											push(@print, $value);
										}
										else {
											print STDERR "Print_Annovar: variant_list: key does not exist. $chr=>$start=>$end=>$ref=>$alt=>$key";
										}
									}
									else {
										print STDERR "Print_Annovar: field_list: key does not exist. $field\n";
									}
								}

								## print total aligner/caller
								foreach my $aligner (sort keys %{$sample_list->{$sample}->{$chr}->{$start}->{$end}->{$ref}->{$alt}}) {
									$aligner_list->{$aligner}++;
									foreach my $caller (sort keys %{$sample_list->{$sample}->{$chr}->{$start}->{$end}->{$ref}->{$alt}->{$aligner}}) {
										$caller_list->{$caller}++;
										$total_score++;
									}
								}
								@print = (@print, $total_score, scalar(keys %{$aligner_list}), scalar(keys %{$caller_list}));

								## print each aligner/caller
								foreach my $aligner (sort keys %{$algorithm_list}) {
									foreach my $caller (sort keys %{$algorithm_list->{$aligner}}) {
										my $flag = 0;
										if(exists $aligner_list->{$aligner} && exists $caller_list->{$caller}) { $flag = 1; }
										push(@print, $flag);
									}
								}

								print OUT join("\t", @print)."\n";
							}
						}
					}
				}
			}
			close(OUT);
			
			print "Output written into $printfile\n";
			
		}
		else {
			print STDERR "Print_Annovar:report_list:$sample:variants:$sample:$key: key does not exist!\n";
		}
	}


	return $variant_list;
}

sub Open_File_Annovar
{
	my ($variant_list, $sample_list, $file, $aligner, $caller) = @_;
	
	open(IN, "<", $file ) or die $!;

	my $i = 0;
	my $pass_flt = 0;
	my $sm_gt_flt;
	my $field_list;
	my $variant_dup_check; ## check whether there are dup variant rows in the same file
	my $sample_index;
	my $format_index = "";
	my $filter_index = "";

	while(my $line = <IN>) {
		chomp($line);
		my @line = split(/\t/, $line);

		if($line =~ m/^Chr\t/) {
			for(my $i=0; $i<@line; $i++) {
				my $field = $line[$i];
				$field =~ s/\s+/_/g;
				# print "$i => $field\n";
				if($field eq "FORMAT") { $format_index = $i; }
				if($field eq "FILTER") { $filter_index = $i; }
				if($format_index ne "" && $i > $format_index) {
					if(exists $sample_index->{$field}) { print STDERR "Open_File_Annovar: sample_index: key already exists. $field\n"; }
					else { $sample_index->{$field} = $i; }
				}

				if(exists $field_list->{$field}) { print STDERR "Open_File_Annovar: field_list: key already exists. $field\n"; }
				else {  $field_list->{$field} = $i; }
			}
		}
		elsif($line =~ m/^chr/) {
			$i++;
			my $filter = $line[$filter_index];

			if($filter ne "PASS") { $pass_flt++; }
			if($pass_flag == 1 && $filter ne "PASS") { next; }

			my $chr = $line[$field_list->{"Chr"}];
			my $start = $line[$field_list->{"Start"}];
			my $end = $line[$field_list->{"End"}];
			my $ref = $line[$field_list->{"Ref"}];
			my $alt = $line[$field_list->{"Alt"}];

			## check dup variant rows in the same file
			if(exists $variant_dup_check->{$chr}->{$start}->{$end}->{$ref}) {
				print STDERR "Open_File_Annovar: variant_dup_check: key already exists! $chr=>$start=>$end=>$ref; file = $file\n";
			}
			else {
				$variant_dup_check->{$chr}->{$start}->{$end}->{$ref} = "";
			}

			## record variant info
			foreach my $key (sort keys %{$field_list}) {
				$variant_list->{$chr}->{$start}->{$end}->{$ref}->{$alt}->{$key} = $line[$field_list->{$key}];
			}

			foreach my $sample (sort keys %{$sample_index}) {
				my $genotype = $line[$sample_index->{$sample}];
				$sm_gt_flt->{$sample} = 0;

				## record only variants in this sample
				## could be "/" or "|" (latter as phased genotype, e.g. freebayes)
				if($genotype =~ m/^(\d+)[\/\|](\d+)/) {
					# print "genotype = $genotype\n";
					if($1 == 0 && $2 == 0) {
						$sm_gt_flt->{$sample}++;
					}
					else {
						## record sample info
						if(exists $sample_list->{$sample}->{$chr}->{$start}->{$end}->{$ref}->{$alt}->{$aligner}->{$caller}) {
							print STDERR "Open_File_Annovar: sample_list: key already exists! $sample=>$chr=>$start=>$end=>$ref=>$alt=>$aligner=>$caller; file = $file\n";
						}
						else {
							$sample_list->{$sample}->{$chr}->{$start}->{$end}->{$ref}->{$alt}->{$aligner}->{$caller} = $file;
						}
					}
				}
				elsif($genotype !~ m/^\./) {
					print STDERR "Open_File_Annovar: genotype not in format. $sample: $genotype; file = $file\n";
				}
			}

		}
		elsif($line =~ m/\S+/) {
			print STDERR "Open_File_Annovar: line not in format. line = $line; file = $file\n"
		}

	}

	close(IN);
	print "Variant total = $i\nVariant PASS = ".($i-$pass_flt)." [ PASS filter is ";
	if($pass_flag) { print "ON"; } else { print "OFF"; }
	print " ]\n";
	print "Variant in each sample: \n##Sample Total PASS\n";
	foreach my $sm (sort keys %{$sm_gt_flt}) {
		print "$sm ".($i-$sm_gt_flt->{$sm})." ";
		if($pass_flag) { print ($i-$pass_flt-$sm_gt_flt->{$sm}); } else { print "NA"; }
		print "\n";
	}
	
	return $variant_list, $sample_list;
}

sub Annovar_Field {
	my $field_list = {
		"Chr" => "Chr",
		"Start" => "Start",
		"End" => "End",
		"Ref" => "Ref",
		"Alt" => "Alt",
		"Func" => "Func.refGene",
		"Gene" => "Gene.refGene",
		"ExonicFunc" => "ExonicFunc.refGene",
		"AAChange" => "AAChange.refGene",
		"dbSNP" => "snp138",
		"ESP6500" => "esp6500siv2_all",
		"Genome1000" => "1000g2014sep_all",
		"ExAC65000" => "exac03",
		"COSMIC" => "cosmic70",
		"CADD" => "CADD_phred",
		"ClinVar" => "clinvar_20140929"
	};

	return $field_list;
}

sub Annovar_Protocol
{
	my $prototol_list = {
		"refGene" => "g",
		"snp138" => "f",
		"esp6500si_all" => "f",
		"1000g2014sep_all" => "f",
		"cosmic70" => "f",
		"clinvar_20140929" => "f",
		"ljb26_all" => "f",
	};

	return $prototol_list;
}

sub Parse_ReportOut
{
	my ($out_list, $report_list) = @_;

	if(exists $out_list->{"data"}) {
		my @data = @{$out_list->{"data"}};
		foreach my $data (@data) {
			my ($sample, $aln_sum, $insert, $qual_by_cycle, $qual_dist, $cov, $variant) = ("") x 7;

			if(exists $data->{"sample"}) { $sample = $data->{"sample"}; }
			if($sample !~ m/\S+/) {
				print STDERR "Parse_ReportOut: out_list:data:sample: key does not exist!\n";
			}
			else {
				if(exists $data->{"fastqc"}->{"inputs"}) { 
					my @inputs = @{$data->{"fastqc"}->{"inputs"}};
					foreach my $input (@inputs) {
						my $rg = "";
						my @keys = qw(paired leftseq rightseq);
						if(exists $input->{"readgroup"}) { $rg = $input->{"readgroup"}; }
						if($rg !~ m/\S+/) {
							print STDERR "Parse_ReportOut: out_list:data:fastqc:input:readgroup: key does not exist!\n";
						}
						else {
							foreach my $key (@keys) {
								$report_list->{$sample}->{"fastqc"}->{$rg}->{$key} = "";
								if(exists $input->{$key}) { 
									$report_list->{$sample}->{"fastqc"}->{$rg}->{$key} = $input->{$key}; 
								}
							}
						}
					}
				}
				if(exists $data->{"alignments"}->{"inputs"}) { 
					my @keys  = qw(alignment_summary_metrics insert_size_metrics quality_by_cycle_metrics quality_distribution_metrics total_coverage);
					foreach my $key (@keys) {
						$report_list->{$sample}->{"alignments"}->{$sample}->{$key} = "";
						if(exists $data->{"alignments"}->{"inputs"}->{$key}) { 
							$report_list->{$sample}->{"alignments"}->{$sample}->{$key} = 
							                                $data->{"alignments"}->{"inputs"}->{$key}; 
						}
					}

				}
				if(exists $data->{"variants"}->{"inputs"}) { 
					my @keys = qw(variant_table);
					foreach my $key (@keys) {
						$report_list->{$sample}->{"variants"}->{$sample}->{$key} = "";
						if(exists $data->{"variants"}->{"inputs"}->{$key}) { 
							$report_list->{$sample}->{"variants"}->{$sample}->{$key} = 
							                    $data->{"variants"}->{"inputs"}->{$key}; 
						}
					}
				}
			}
		}

	}
	else {
		print STDERR "Parse_ReportOut:out_list:data: key does not exist!\n";
	}


	return $report_list;
}

sub Print_Menu
{
	my $PROG = shift;

	my $menu = qq~
# This script prepares and processes fastqc/aln/annovar files as input for pipeline report utility.

Usage: 
 $PROG -i|--in <report_in_yaml> -o|--out <report_out_yaml> -d|--dir <output_directory> [ -p|--pass]

Options:
 [-i|--in]      : YAML file with information of input files generated by the pipeline.
 [-o|--out]     : YAML file with information of output files for the project report.
 [-d|--dir]     : Output directory to print output files.
 [-p|--pass]    : Print PASS variants only.

Example: 
 $PROG -i LCAexome.report.in.yaml -o LCAexome.report.out.yaml -d LCAexome_report -p

Release: 2015-02-27
Contact: Riyue Bao <rbao\@uchicago.edu>

~;

	return $menu;
}

sub Update_Progress
{
	my $progress = shift;
	
	print "[ ",Local_Time(), " ] $progress\n";
	
	return $progress;
}

sub Parse_Inputfiles
{
	my @files = @_;

	my @files_new;
	
	if(@files == 1 && $files[0] =~ m/^(\S*)\*(\S*)$/)
	{
		@files_new = <$1*$2>;
	}
	elsif(@files > 1)
	{
		foreach my $file(@files)
		{
			if($file =~ m/^(\S*)\*(\S*)$/)
			{
				@files_new = (@files_new, <$1*$2>);
			}
			else
			{
				push(@files_new, $file);
			}
		}
	}
	else
	{
		@files_new = @files;
	}

	return @files_new;
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

