
package ExScaliburGMD::PrintIVC;

use strict;
use warnings;
use Exporter qw(import);
use Cwd qw(abs_path);
use File::Basename;
use lib dirname(abs_path($0)).'/lib/perl5';
use ExScaliburGMD::Util qw(Empty_Value);

our @EXPORT = qw(Print_IVC);
our @EXPORT_OK = qw(Print_IVC);

## ---------------------------------------

sub Print_IVC
{
  my ($config_list, $output_dir, $output_file) = @_;

  open(OUT, ">", "$output_dir/$output_file") or die $!;
  
  print OUT "\n[user]\n";
  if(exists $config_list->{"pipeline"}->{"ivc"}->{"config"}) {
    foreach my $key (sort keys %{$config_list->{"pipeline"}->{"ivc"}->{"config"}}) {
      my $value = "";
      $value = Empty_Value($value, $config_list->{"pipeline"}->{"ivc"}->{"config"}->{$key});
      print OUT "$key = $value\n";
    }
  }
  else {
    print STDERR "Print_IVC: config_list: key does not exist. pipeline=>ivc=>config\n";  
  }

  close(OUT);

  return $config_list;

}


## ---------------------------------------

1;

