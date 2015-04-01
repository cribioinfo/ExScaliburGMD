#!/usr/bin/perl -w

=head1 LICENSE

Collect_results.pl

Copyright (C) 2013-2015 Center for Research Informatics, The University of Chicago

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

############################################################################################
# Riyue Bao 03/10/2015
# This script collects result files to archive directory and prepare output.tar.gz for downloading.
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
use File::Copy qw(copy move);

## options
use vars qw($opt_tar);
use vars qw($opt_reverse);

##############################
# Main
##############################

## version
my $version = "1.0.0";
my $release_date = "2015-03-10";
my $desc = "Collect project config and result files into archive directory";

## initialize 
my $config_list;
my $program;
my $menu;
my $command;
my $file = "";
my $output = "myProject";
my $output_file = "";
my $project = "";
my $archive_dir = "";
my $tar_flag = 0;
my $reverse_flag = 0;
my @dirs = qw(alignment alignment_metrics  configs  project_report  qc_reports  variant_annotation  variant_calls);

## print menu
$program = basename($0);
$menu = Print_Menu($program, $version, $release_date, $desc);
if(@ARGV == 0) { print "\n$menu\n"; exit; }

## print command and info
$command = Get_Opt($command);
$command = Print_Info($command);

## run start
$output_file = "$output.slim.tar.gz";

Update_Progress("Read pipeline archive configuration from $file");
$config_list = YAML::Tiny::LoadFile($file);

# print Dumper($config_list);
# print $config_list->{"output"}."\n";
# YAML::Tiny::DumpFile("test.out",$config_list);

if(exists $config_list->{"project"}) { 
	$project = $config_list->{"project"}; 
	print "project = $project\n" ;
}
else { 
	print STDERR "Compress_Files: config_list: key does not exist. key = project\n"; 
}

if(exists $config_list->{"directory"}) { 
	$archive_dir = $config_list->{"directory"}; 
	print "archive directory = $archive_dir\n" ;
}
else { 
	print STDERR "Compress_Files: config_list: key does not exist. key = directory\n"; 
}

## copy/move files
$config_list = Archive_Files($config_list);

## compress archive dir
if($tar_flag) {
	Update_Progress("Compressing \"$archive_dir\" into $output_file");
	print "Note: Alignment BAM and BED files will be excluded!\n";
	$config_list = Compress_Files($config_list, $output_file, $archive_dir, \@dirs);
}

## run end
Update_Progress("Program run finished\n");

##############################
# Functions
##############################

sub Archive_Files
{
	my $config_list = shift;

	my $file_count;

	foreach my $type (qw(file dir)) {
		foreach my $action (qw(copy move)) {
			$file_count->{$type}->{$action} = 0;
		}
	}

	if(exists $config_list->{"data"}) {
		my @data = @{$config_list->{"data"}};

		foreach my $data (@data) {
			foreach my $dir (sort keys %{$data}) {
				if($dir eq "sample") {
					print "sample = ".$data->{"sample"}."\n";
				}
				else {
					my @info = @{$data->{$dir}};
					foreach my $info (@info) {
						## keys  = to from action type
						my $to = $info->{"to"};
						my $from = $info->{"from"};
						my $action = $info->{"action"};
						my $type = $info->{"type"};
						
						$file_count->{$type}->{$action}++;

						## start moving/copying files...
						$from = Action_File($from, $to, $action, $type);

						if($reverse_flag) {
							## move files back if --reverse is on
							Update_Progress("Reverse collection procedure: moving files back to their original location");
							$to = Reverse_File($from, $to, $action, $type);
						}
					}
				}
			}
		}

	}
	else {
		print STDERR "Archive_Files: config_list: key does not exist. key = directory\n";
	}

	return $config_list;
}

sub Reverse_File
{
	my ($from, $to, $action, $type) = @_;

	if($type eq "file") {
		if(-f $to) {
			if($action eq "move") {
				# `mv $to $from`;
			}
		}	
	}
	elsif($type eq "dir") { 
		if(-d $to) {
			if($action eq "move") {
				# `mv $to/* $from/`;
			}
		}
	}

	return $to;
}

sub Action_File
{
	my ($from, $to, $action, $type) = @_;

	if($type eq "file") {
		if(-f $from) {
			if($action eq "copy") {
				# `cp -p $from $to`;
			}
			elsif($action eq "move") {
				# `mv $from $to`;
			}
			else {
				print STDERR "Archive_Files: action not in format. Must be move or copy. action = $action\n";
			}
		}
		else {
			print STDERR "Archive_Files: file does not exist. file = $from\n";
		}	
	}
	elsif($type eq "dir") { ## process all files in the original directory (NOT the directory itself!)
		if(-d $from) {
			if($action eq "copy") {
				# `cp -p $from/* $to/`;
			}
			elsif($action eq "move") {
				# `mv $from/* $to/`;
			}
			else {
				print STDERR "Archive_Files: action not in format. Must be move or copy. action = $action\n";
			}
		}
		else {
			print STDERR "Archive_Files: dir does not exist. dir = $from\n";
		}	
	}
	else {
		print STDERR "Archive_Files: type not in format. Must be file or dir. type = $type\n"; 
	}

	return $from;
}

sub Compress_Files
{
	my ($config_list, $output_file, $archive_dir, $dirs_ref) = @_;

	my $dir_string = "";
	if($archive_dir ne "") {
		foreach my $dir (@{$dirs_ref}) { 
			if($dir ne "alignment") { $dir_string .= " $archive_dir/$dir"; }
		}
		`tar zcvf $output_file $dir_string`;
	}
	else {
		print STDERR "Compress_Files: Archive directory is missing. Check key \"directory\" in $file: Is it missing or has empty value? Tar compression failed!\n";
	}
	
	return $config_list;
}

sub Get_Opt
{
	my $command = shift;

	GetOptions( 
		"file|f:s" => \$file,
		"output|o:s" => \$output,
		"tar|t",
		"reverse|r"
	) or die $!;

	## check
	if($file eq "") { print STDERR "Input file missing! \n$menu\n"; exit(1); }	
	if($opt_tar) { $tar_flag = 1; }
	if($opt_reverse) { $reverse_flag = 1;}

	## command
	$command = "$program --file $file --output $output";
	if($tar_flag) { $command .= " --tar"; }
	if($reverse_flag) { $command .= " --reverse"; }

	return $command;
}

sub Print_Info
{
	my $command = shift;

	print "\n[COMMANDS]\n$command\n\n[PROGRESS]\n";

	return $command;
}

sub Print_Menu
{
	my $program = shift;

	my $menu = qq~
Program: $program
Description: $desc
Version: $version (release $release_date)
Contact: Riyue Bao <rbao\@uchicago.edu>

Usage: $program -f|--file <yaml> [OPTIONS]

Required:
 [-f|--file]      : YAML file with information of file names and paths to be collected.
                    Default is alignment BAM and variant VCF and annotation files will be 
                    moved due to their large size, and other files will be copied. 
                    You may specify which files to be copied or moved in the YAML file. 
             Note : Upon "copy", the original copy of the file will stay intact.
             Note : Upon "move", the file will be removed from the original location.

Optional: 
 [-t|--tar]       : Generate a tar.gz file from the output archive directory (alignment
                    files excluded). Default is off.
 [-o|--output]    : Output file prefix. Default is "myProject". Will be use to name the 
                    tar.gz file "myProject.slim.tar.gz" if --tar is on.
             Note : "slim" indicates the tar.gz file does not include alignment files.                     
 [-r|--reverse]   : Reverse the collection procedure and move all files back to their 
                    original location. Use this option with extra caution!
             Note : The original directory where the file comes from must stay intact, 
                    as well as all parental directories.                    

Example: 
 $program -f LCAexome.archive.yaml -t -o LCAexomeProj/archive/LCAexome_archive

Note:
 The script does not generate new directories or files. All necessary read and write directories 
 must be present prior to running the code.
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

