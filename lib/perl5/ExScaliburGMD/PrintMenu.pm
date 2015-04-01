
package ExScaliburGMD::PrintMenu;

use strict;
use warnings;
use Exporter qw(import);

our @EXPORT = qw(Print_Menu);
our @EXPORT_OK = qw(Print_Menu);

## ---------------------------------------

sub Print_Menu
{
  my ($program, $author, $email, $version, $release_date, $license, $desc) = @_;

  my $menu = qq~
Program     : $program
Description : $desc
Version     : $version (release $release_date)
License     : $license
Contact     : $author <$email>

Usage: $program -m|--metadata <metadata_table> -c|--config <pipeline_config> [OPTIONS]

Required:
 [-m|--metadata]    : Sample metadata table in tab-delimited format.
 [-c|--config]      : Pipeline configuration file in YAML format.

Optional:
 [-p|--project]     : Project title. (default: myProject)
 [-d|--projdir]     : Project output directory. (default: the current directory)
               Note : Must be full absolute path!
 [-t|--threads]     : Number of threads. (default: 1)
 [-q|--mapq]        : Filter out reads with MapQ (mapping quality) lower than this 
                      threshold in the alignment. (default: 0)
 [--split]          : Activate scatter-gather mode. (default: off) [*].
 [--force]          : Include all samples and readgroup regardless of whether it 
                      fails QC. (default: off)
 [--tree]           : tree executable (default: tree). Must be in system path.
               Note : Tree is a recursive directory listing program.
                      If specified, project directory structure will be printed to 
                      project.tree.txt.
Developer: 
 [--local]          : Running on a single computer.
 [--ssh]            : A server farm or a bunch of desktops or servers without a
                      a workload management system.
 [--cluster]        : Running on a cluster (HPC or cloud).
               Note : If no platform options are specified, the pipeline will run -local.
 [--moab]           : Running on a MOAB/PBS cluster.
 [--pbs]            : Running on a PBS cluster.
 [--sge]            : Running on a SGE cluster.
 [-y|--retry]       : Number of times to retry a failing tasks. (default: 0; maximum: 10)
 [-l|--log]         : Log all tasks (do not delete tmp files). (default: off)
 [--bdscfg]         : BigDataScript config file. Users may provide a custom BDS config.
                      If specified, the pipeline will copy this file to BDS config in
                      the project directory.
                      If not specified, the pipeline will generate BDS config from 
                      default settings.
Notes: 
 [*] scatter-gather : If activated, target regions will be split into even-sized bins 
                      and variants will be called from each bin. At the end, results of
                      each bin are merged together to generate one complete set of 
                      variants for the whole exome. The [split] mode speeds up the analysis
                      considerably provided with sufficient computational resources.

Example: 
 $program -m LCAexome.metadata.txt -c LCAexome.pipeline.yaml -p LCAexome -d pathTo/LCAexomeProj -t 4 -q 30 --force --split --tree /group/bioinformatics/software/tree-1.6.0/tree -cluster -sge -y 3
~;

  return $menu;
}

## ---------------------------------------

1;
