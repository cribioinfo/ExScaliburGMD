
package ExScaliburGMD::PrintReadgroup;

use strict;
use warnings;
use Exporter qw(import);
use Cwd qw(abs_path);
use File::Basename;
use lib dirname(abs_path($0)).'/lib/perl5';
use ExScaliburGMD::Util qw(Empty_Value Value_And_Exit Value_And_Ctd);

our @EXPORT = qw(Print_Readgroup);
our @EXPORT_OK = qw(Print_Readgroup);

## ---------------------------------------

sub Print_Readgroup
{
  my ($sample_list, $config_list, $proj_dir, $project, $sample, $library, $readgroup) = @_;

  ## initialize
  my $output_file = "$sample.$readgroup.cfg";
  my $output_dir = "$proj_dir/configs/$project\_samples/$sample";
  my $RGstring = "";
  my $genome = "";
  my $flavor = "";
  my $location = "";
  my $read_length = 0;
  my $paired_flag = 0;
  my @print_files;
  my @keys = qw(Seqfile1 Seqfile2);

  $RGstring = Value_And_Exit($RGstring, $sample_list->{$sample}->{$library}->{$readgroup}, "RGstring");
  $genome = Value_And_Exit($RGstring, $sample_list->{$sample}->{$library}->{$readgroup}, "Genome");
  $flavor = Value_And_Exit($RGstring, $sample_list->{$sample}->{$library}->{$readgroup}, "Flavor");
  $location = Value_And_Ctd($RGstring, $sample_list->{$sample}->{$library}->{$readgroup}, "Location");

  if($flavor =~ m/^(\d+)x(\d+)$/) {
    $read_length = $2;
    if($1 == 2) { $paired_flag = 1; }
  }
  if($read_length == 0) {
    print STDERR "Print_Readgroup: sample_list: Flavor is not in format (1 or 2 x readlength). Currently it is \"$flavor\". Program terminated!\n";
    exit;
  }

  ## print attributes
  open(OUT, ">", "$output_dir/$output_file") or die $!;

  ## -----
  print OUT "\n## project/sample/readgroup information\n";
  print OUT "project = $project\n";
  print OUT "sample = $sample\n";
  print OUT "library = $library\n";
  print OUT "readgroup = $readgroup\n";
  print OUT "rg_string = \"$RGstring\"\n";
  print OUT "genome_assembly = $genome\n";
  print OUT "read_length = $read_length\n";
  print OUT "paired = $paired_flag\n";

  ## -----
  @print_files = ();
  print OUT "\n## seqfiles\n";
  foreach my $key (sort @keys) {
    if(exists $sample_list->{$sample}->{$library}->{$readgroup}->{$key}) {
      push(@print_files, "$location/".$sample_list->{$sample}->{$library}->{$readgroup}->{$key});
    }
  }
  print OUT "seqfiles = ",join(",", @print_files)."\n";

  ## -----
  @print_files = ();
  print OUT "\n## qc reports\n";
  foreach my $key (sort @keys) {
    if(exists $sample_list->{$sample}->{$library}->{$readgroup}->{$key}) {
      my $seqfile = $sample_list->{$sample}->{$library}->{$readgroup}->{$key};
      $seqfile =~ s/.fastq.gz$//;
      $seqfile =~ s/.fq.gz$/.fq/;
      push(@print_files, "$proj_dir/results/$project\_samples/$sample/qc_reports/$seqfile\_fastqc.zip");
    }
  }
  print OUT "fastqc_files = ",join(",", @print_files)."\n";

  ## -----
  @print_files = ();
  ## paired-end
  if($paired_flag == 1) {
    print OUT "\n## seqprep files\n";
    foreach my $key (qw(R1 R2 unpaired)) {
      push(@print_files, "$proj_dir/results/$project\_samples/$sample/clean_reads/$readgroup.$key.clip.fq.gz");
    }
    print OUT "seqprep_files = ",join(",", @print_files)."\n";
  }
  ## single-end
  else {
    print OUT "\n## cutadapt files\n";
    print OUT "cutadapt_files = $proj_dir/results/$project\_samples/$sample/clean_reads/$readgroup.clip.fq.gz\n";
  }
  
  close(OUT);

  return $sample_list;
}

## ---------------------------------------

1;

