
package ExScaliburGMD::BuildReportYaml;

use strict;
use warnings;
use Exporter qw(import);
use Cwd qw(abs_path);
use File::Basename;
use YAML::Tiny;
use Data::Dumper;
use lib dirname(abs_path($0)).'/lib/perl5';
use ExScaliburGMD::Util qw(Empty_Value Value_And_Exit Value_And_Ctd Add_Key BAM_Suffix Set_Fastqc_File);

our @EXPORT = qw(Build_Report_Yaml);
our @EXPORT_OK = qw(Build_Report_Yaml);

## ---------------------------------------

sub Build_Report_Yaml
{
  my ($sample_list, $config_list, $proj_dir, $project) = @_;

  ## initialize
  my $aligner_list;
  my $caller_list;
  my $suffix = "";
  my $fastqc_version = undef;
  my $annovar_version = undef;
  my @keys = qw(alignment_summary_metrics insert_size_metrics quality_by_cycle_metrics quality_distribution_metrics total_coverage);

  ## decide refined bam filename based on pipeline module flags
  $suffix = BAM_Suffix($suffix, $config_list->{"pipeline"}->{"flags"}->{"modules"});

  ## retrieve aligner/caller list
  $aligner_list = Add_Key($config_list->{"pipeline"}->{"flags"}, "aligners", $aligner_list);
  $caller_list = Add_Key($config_list->{"pipeline"}->{"flags"}, "callers", $caller_list);

  ## print yaml
  $config_list = Print_Report_Yaml_In($sample_list, $config_list, $proj_dir, $project, $aligner_list, $caller_list, $suffix, $fastqc_version, $annovar_version, \@keys);

  $config_list = Print_Report_Yaml_Out($sample_list, $config_list, $proj_dir, $project, $aligner_list, $caller_list, $suffix, $fastqc_version, $annovar_version, \@keys);

  return $config_list;
}

sub Print_Report_Yaml_In
{
  my ($sample_list, $config_list, $proj_dir, $project, $aligner_list, $caller_list, $suffix, $fastqc_version, $annovar_version, $keys_ref) = @_;

  my $report_list;
  my $key_list;
  my $output_file = "$proj_dir/configs/$project\_project/$project.report.in.yaml";
  my $paired_flag = 0;
  my @data;
  my @variants;
  my @keys = @{$keys_ref};

  ## set hash for printing
  foreach my $sample (sort keys %{$sample_list}) {
    my @rgs;
    my $paired_flag_sm = 0;

    ## build readgroup list
    foreach my $library (sort keys %{$sample_list->{$sample}}) {

      ## record rgs
      foreach my $readgroup (sort keys %{$sample_list->{$sample}->{$library}}) {
        my $leftseq = undef;
        my $rightseq = undef;
        my $paired_flag = 0;
        my $flavor = "";

        $leftseq = Set_Fastqc_File($leftseq, $sample_list->{$sample}->{$library}->{$readgroup}->{"Seqfile1"});
        $rightseq = Set_Fastqc_File($rightseq, $sample_list->{$sample}->{$library}->{$readgroup}->{"Seqfile2"});
        $flavor = Empty_Value($flavor, $sample_list->{$sample}->{$library}->{$readgroup}->{"Flavor"});
        
        if($flavor =~ m/^2x\d+$/) { $paired_flag = 1; }
        elsif($flavor =~ m/^1x\d+$/) { $paired_flag = 0; }
        else {
          print STDERR "Build_Report_Yaml:In: sample_list: Flavor is not in format (1 or 2 x readlength). Currently it is \"$flavor\". Program terminated!\n";
          exit;
        } 

        ## if at least one of the rgs is PE, the sample is set to PE
        if($paired_flag) { $paired_flag_sm = 1; }
        
        ## build hash
        my $rg_element = {
          "path" => "$proj_dir/results/$project\_samples/$sample/qc_reports",
          "readgroup" => $readgroup,
          "leftseq" => $leftseq,
          "rightseq" => $rightseq,
          "paired" => $paired_flag
        };
        
        push(@rgs, $rg_element);
      }
    }

    ## set alignment attributes
    foreach my $key (@keys) {
      my @files;
      foreach my $aligner (sort keys %{$aligner_list}) {
        push(@files, "$sample.$aligner.merged".$suffix.".bam.$key");
      }
      $key_list->{$key} = join(",", @files);

      ## single-end reads will not have picard insert_size output file
      if($key eq "insert_size_metrics" && $paired_flag_sm == 0) {
        $key_list->{$key} = undef;
      }
    }

    ## build bash
    my $data_element = {
      "fastqc" => {
        "version" => $fastqc_version, 
        "inputs" => [@rgs]
      },
      "alignments" => {
        "inputs" => {
          "path" => "$proj_dir/results/$project\_samples/$sample/alignment_metrics",
          "alignment_summary_metrics" => $key_list->{"alignment_summary_metrics"},
          "insert_size_metrics" => $key_list->{"insert_size_metrics"},
          "quality_by_cycle_metrics" => $key_list->{"quality_by_cycle_metrics"},
          "quality_distribution_metrics" => $key_list->{"quality_distribution_metrics"},
          "total_coverage" => $key_list->{"total_coverage"}
        }
      },
      "sample" => $sample
    };
    
    push(@data, $data_element);
  } 

  ## build hash
  foreach my $aligner (sort keys %{$aligner_list}) {
    foreach my $caller (sort keys %{$caller_list}) {
      my $variant_element = {
        "path" => "$proj_dir/results/$project\_multisample/variant_annotation",
        "variant_table" => "$project.$aligner.$caller.flt.anno.exonic.txt"
      };
      push(@variants, $variant_element);
    }
  }

  $report_list = {
    "project" => $project,
    "data" => [@data],
    "multi" => {
      "variants" => {
        "inputs" => [@variants],
        "version" => $annovar_version
      }
    }
  };

  # print Dump($report_list), "\n";
  YAML::Tiny::DumpFile($output_file, $report_list);

  return $config_list;
}

sub Print_Report_Yaml_Out
{
  my ($sample_list, $config_list, $proj_dir, $project, $aligner_list, $caller_list, $suffix, $fastqc_version, $annovar_version, $keys_ref) = @_;

  my $report_list;
  my $key_list;
  my $output_file = "$proj_dir/configs/$project\_project/$project.report.out.yaml";
  my $paired_flag = 0;
  my @data;
  my @variants;
  my @keys = @{$keys_ref};

  ## set hash for printing
  foreach my $sample (sort keys %{$sample_list}) {
    my @rgs;
    my $paired_flag_sm = 0;

    ## build readgroup list
    foreach my $library (sort keys %{$sample_list->{$sample}}) {

      ## record rgs
      foreach my $readgroup (sort keys %{$sample_list->{$sample}->{$library}}) {
        my $leftseq = undef;
        my $rightseq = undef;
        my $paired_flag = 0;
        my $flavor = "";

        $leftseq = Set_Fastqc_File($leftseq, $sample_list->{$sample}->{$library}->{$readgroup}->{"Seqfile1"});
        $rightseq = Set_Fastqc_File($rightseq, $sample_list->{$sample}->{$library}->{$readgroup}->{"Seqfile2"});
        $flavor = Empty_Value($flavor, $sample_list->{$sample}->{$library}->{$readgroup}->{"Flavor"});
        
        if($flavor =~ m/^2x\d+$/) { $paired_flag = 1; }
        elsif($flavor =~ m/^1x\d+$/) { $paired_flag = 0; }
        else {
          print STDERR "Build_Report_Yaml:Out: sample_list: Flavor is not in format (1 or 2 x readlength). Currently it is \"$flavor\". Program terminated!\n";
          exit;
        } 

        ## single-end does not have rightseq (should stay as undef)
        if(defined $leftseq) { $leftseq = "fastqc/$leftseq"; }
        if(defined $rightseq) { $rightseq = "fastqc/$rightseq"; }

        ## if at least one of the rgs is PE, the sample is set to PE
        if($paired_flag) { $paired_flag_sm = 1; }
        
        ## build hash
        my $rg_element = {
          "readgroup" => $readgroup,
          "leftseq" => $leftseq,
          "rightseq" => $rightseq,
          "paired" => $paired_flag
        };
        
        push(@rgs, $rg_element);
      }
    }

    ## set alignment attributes
    ## all aligners in one file
    foreach my $key (@keys) {
      $key_list->{$key} = "aln/$sample.$key.txt";

      ## single-end reads will not have picard insert_size output file
      if($key eq "insert_size_metrics" && $paired_flag_sm == 0) {
        $key_list->{$key} = undef;
      }
    }

    ## build bash
    my $data_element = {
      "fastqc" => {
        "version" => $fastqc_version, 
        "inputs" => [@rgs]
      },
      "alignments" => {
        "inputs" => {
          "alignment_summary_metrics" => $key_list->{"alignment_summary_metrics"},
          "insert_size_metrics" => $key_list->{"insert_size_metrics"},
          "quality_by_cycle_metrics" => $key_list->{"quality_by_cycle_metrics"},
          "quality_distribution_metrics" => $key_list->{"quality_distribution_metrics"},
          "total_coverage" => $key_list->{"total_coverage"}
        }
      },
      "variants" => {
        "inputs" => {
          "variant_table" => "annovar/$sample.variants.tsv"
        }
      },
      "sample" => $sample
    };
    
    push(@data, $data_element);
  } 

  $report_list = {
    "project" => $project,
    "data" => [@data],
  };

  # print Dump($report_list), "\n";
  YAML::Tiny::DumpFile($output_file, $report_list);

  return $config_list;
}

## ---------------------------------------

1;

