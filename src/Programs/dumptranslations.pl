#!/usr/local/bin/perl


use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my $host   = 'kaka.sanger.ac.uk';
my $user   = 'anonymous';
my $dbname = 'homo_sapiens_core_12_31';
my @chromosomes;
my $file = 'stdout';

my $path = 'NCBI31';
my $port = 3306;

$| = 1;

&GetOptions(
  'host:s'   => \$host,
  'user:s'   => \$user,
  'dbname:s' => \$dbname,
  'path:s'   => \$path,
  'port:n'   => \$port,
  'chromosomes:s' => \@chromosomes,
  'file:s' => \$file,
);

if (scalar(@chromosomes)) {
  @chromosomes = split(/,/,join(',',@chromosomes));
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $host,
  -user   => $user,
  -port   => $port,
  -dbname => $dbname
);
$db->assembly_type($path);

my $sa   = $db->get_SliceAdaptor();

my $chrhash = get_chrlengths($db);

#filter to specified chromosome names only
if (scalar(@chromosomes)) {
  foreach my $chr (@chromosomes) {
    my $found = 0;
    foreach my $chr_from_hash (keys %$chrhash) {
      if ($chr_from_hash =~ /^${chr}$/) {
        $found = 1;
        last;
      }
    }
    if (!$found) {
      print "Didn't find chromosome named $chr in database $dbname\n";
    }
  }
  HASH: foreach my $chr_from_hash (keys %$chrhash) {
    foreach my $chr (@chromosomes) {
      if ($chr_from_hash =~ /^${chr}$/) {next HASH;}
    }
    delete($chrhash->{$chr_from_hash});
  }
}


$| = 1;
if ($file ne "stdout") {
  open FP,">$file";
} else {
  open FP,">-";
}

foreach my $chr (reverse sort bychrnum keys %$chrhash) {
  my $chrstart = 1;
  my $chrend   = $chrhash->{$chr};

  $chrend=230000000;
  print STDERR "Chr $chr from 1 to " . $chrend. "\n";

  my $slice = $sa->fetch_by_chr_start_end($chr,$chrstart,$chrend);

  print "Fetching genes\n";
  my $genes = $slice->get_all_Genes;

  print "Done fetching genes (fetched " . scalar(@$genes) .")\n";

# Go through the genes sorting the transcripts
  print "Starting sorting transcripts\n";
  foreach my $gene (@$genes) {
    foreach my $trans (@{$gene->get_all_Transcripts}) {
      $trans->sort;
    }
  }
  print "Done sorting transcripts\n";
  print "Starting translating transcripts\n";
  foreach my $gene (@$genes) {
    print "Gene " . $gene->stable_id . "\n";
    foreach my $trans (@{$gene->get_all_Transcripts}) {
      print " Transcript " . $trans->stable_id . "\n";
      my $pep = $trans->translate;
      print $pep->seq . "\n";
    }
  }
  print "Done translating transcripts\n";
}
print "Done\n";

sub bychrnum {

  my @awords = split /_/, $a;
  my @bwords = split /_/, $b;

  my $anum = $awords[0];
  my $bnum = $bwords[0];

  #  if ($anum !~ /^chr/ || $bnum !~ /^chr/) {
  #    die "Chr name doesn't begin with chr for $a or $b";
  #  }

  $anum =~ s/chr//;
  $bnum =~ s/chr//;

  if ($anum !~ /^[0-9]*$/) {
    if ($bnum !~ /^[0-9]*$/) {
      return $anum cmp $bnum;
    } else {
      return 1;
    }
  }
  if ($bnum !~ /^[0-9]*$/) {
    return -1;
  }

  if ($anum <=> $bnum) {
    return $anum <=> $bnum;
  } else {
    if ($#awords == 0) {
      return -1;
    } elsif ($#bwords == 0) {
      return 1;
    } else {
      return $awords[1] cmp $bwords[1];
    }
  }
}

sub get_chrlengths {
  my $db   = shift;
  my $type = shift;

  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "get_chrlengths should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }

  my $ca = $db->get_ChromosomeAdaptor();
  my $chroms = $ca->fetch_all();

  my %chrhash;

  foreach my $chr (@$chroms) {
    $chrhash{$chr->chr_name} = $chr->length;
  }
  return \%chrhash;
}


