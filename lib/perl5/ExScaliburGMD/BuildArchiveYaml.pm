
package ExScaliburGMD::BuildArchiveYaml;

use strict;
use warnings;
use Exporter qw(import);
use Cwd qw(abs_path);
use File::Basename;
use YAML::Tiny;
use Data::Dumper;
use lib dirname(abs_path($0)).'/lib/perl5';
use ExScaliburGMD::Util qw(Empty_Value Value_And_Exit Value_And_Ctd Add_Key BAM_Suffix Set_Fastqc_File);

our @EXPORT = qw(Build_Archive_Yaml);
our @EXPORT_OK = qw(Build_Archive_Yaml);

## ---------------------------------------

sub Build_Archive_Yaml
{
  my ($sample_list, $config_list, $proj_dir, $project) = @_;

  ## initialize
  my $report_list;
  my $aligner_list;
  my $caller_list;
  my $suffix = "";
  my $output_dir = "$proj_dir/archive/$project\_archive";
  my $output_file = "$proj_dir/configs/$project\_project/LCAexome.archive.yaml";
  my @data;

  ## decide refined bam filename based on pipeline module flags
  $suffix = BAM_Suffix($suffix, $config_list->{"pipeline"}->{"flags"}->{"modules"});

  ## retrieve aligner/caller list
  $aligner_list = Add_Key($config_list->{"pipeline"}->{"flags"}, "aligners", $aligner_list);
  $caller_list = Add_Key($config_list->{"pipeline"}->{"flags"}, "callers", $caller_list);

  ## set hash for printing
  foreach my $sample (sort keys %{$sample_list}) {

    ## initialize
    my @cfgs;
    my @alns;
    my @metrics;
    my @qc;
    my @variants;
    my @annotation;
    my $fn;

    ## assign values
    ## -----
    my $cfg = {
      "action" => "copy",
      "type" => "dir",
      "from" => "$proj_dir/configs",
      "to" => "$output_dir/configs"
    };
    push(@cfgs, $cfg);

    ## -----
    foreach my $aligner (sort keys %{$aligner_list}) {
      $fn = "$sample.$aligner.merged$suffix.loci.callable.bed";
      my $aln = {
        "action" => "move",
        "type" => "file",
        "from" => "$proj_dir/$project\_samples/$sample/alignment/$fn",
        "to" => "$output_dir/alignment/$fn"
      };
      push(@alns, $aln);
    };

    ## -----
    my $metrics = {
      "action" => "copy",
      "type" => "dir",
      "from" => "$proj_dir/results/$project\_samples/$sample/alignment_metrics",
      "to" => "$output_dir/alignment_metrics"
    };
    push(@metrics, $metrics);

    ## -----
    my $qc = {
      "action" => "copy",
      "type" => "dir",
      "from" => "$proj_dir/results/$project\_samples/$sample/qc_reports",
      "to" => "$output_dir/qc_reports"
    };
    push(@qc, $qc);
    foreach my $string (qw(qc_summary qc_flag)) {
      $fn = "$sample.$string.tsv";
      $qc = {
        "action" => "copy",
        "type" => "file",
        "from" => "$proj_dir/results/$project\_samples/$sample/$fn",
        "to" => "$output_dir/qc_reports/$fn"
      };
      push(@qc, $qc);
    };

    ## -----
    foreach my $aligner (sort keys %{$aligner_list}) {
      foreach my $caller (sort keys %{$caller_list}) {  
        foreach my $string (qw(vcf vcf.idx)) {
          $fn = "$project.$aligner.$caller.flt.$string";
          my $var = {
            "action" => "move",
            "type" => "file",
            "from" => "$proj_dir/results/$project\_multisample/variant_calls/$fn",
            "to" => "$output_dir/variant_calls/$fn"
          };
          push(@variants, $var);
        }
      }
    };

    ## -----
    foreach my $aligner (sort keys %{$aligner_list}) {
      foreach my $caller (sort keys %{$caller_list}) {
        foreach my $string (qw(anno anno.exonic)) {
          $fn = "$project.$aligner.$caller.flt.$string.txt";
          my $anno = {
          "action" => "move",
          "type" => "file",
          "from" => "$proj_dir/results/$project\_multisample/variant_calls/$fn",
          "to" => "$output_dir/variant_calls/$fn"
          };
          push(@annotation, $anno);
        }
      }
    };

    ## build bash
    my $data_element = {
      "configs" => [@cfgs],
      "alignments" => [@alns],
      "alignment_metrics" => [@metrics],
      "qc_reports" => [@qc],
      "variant_calls" => [@variants],
      "variant_annotation" => [@annotation],
      "sample" => $sample
    };
    
    push(@data, $data_element);
  } 

  $report_list = {
    "project" => $project,
    "data" => [@data],
    "directory" => $output_dir
  };

  # print Dump($report_list), "\n";
  YAML::Tiny::DumpFile($output_file, $report_list);


  return $config_list;
}



## ---------------------------------------

1;

