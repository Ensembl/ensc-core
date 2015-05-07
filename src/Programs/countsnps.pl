#!/usr/local/bin/perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::Tools::CodonTable;
use Getopt::Long;

my $host   = 'ecs2b';
my $user   = 'ensro';
my $dbname = 'steve_chr6_new_typed_snp';
my @chromosomes;
my $file = 'stdout';
my @snptypes;
my @defgenetypes = ('Known','Novel_CDS');
my @genetypes;

my $path = 'SANGER_13';
my $port = 3306;

$| = 1;

&GetOptions(
  'host:s'   => \$host,
  'user:s'   => \$user,
  'dbname:s' => \$dbname,
  'path:s'   => \$path,
  'port:n'   => \$port,
  'chromosomes:s' => \@chromosomes,
  'snptypes:s' => \@snptypes,
  'genetypes:s' => \@genetypes,
  'file:s' => \$file,
);

if (scalar(@chromosomes)) {
  @chromosomes = split(/,/,join(',',@chromosomes));
}

if (scalar(@snptypes)) {
  @snptypes = split(/,/,join(',',@snptypes));
}

if (scalar(@genetypes)) {
  @genetypes = split(/,/,join(',',@genetypes));
} else {
  @genetypes = @defgenetypes;
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $host,
  -user   => $user,
  -port   => $port,
  -dbname => $dbname
);
$db->assembly_type($path);

my $sa   = $db->get_SliceAdaptor();
my $ga   = $db->get_GeneAdaptor();
my $dafa = $db->get_DnaAlignFeatureAdaptor();

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
  print STDERR "Chr $chr from 1 to " . $chrhash->{$chr}. "\n";
  my $chrstart = 1;
  my $chrend   = $chrhash->{$chr};

  my $slice = $sa->fetch_by_chr_start_end($chr,$chrstart,$chrend);

  print "Fetching genes\n";
  my %genes_hash;
  my @genes;
  foreach my $genetype (@genetypes) {
    $genes_hash{$genetype} = $slice->get_all_Genes_by_type($genetype);
    print "Got " . scalar(@{$genes_hash{$genetype}}) . " $genetype genes\n";
    push @genes,@{$genes_hash{$genetype}};
  }

  print "Done fetching genes (fetched " . scalar(@genes) .")\n";

  @genes = sort {$a->start <=> $b->start} @genes;


  print "Fetching SNPs\n";
  my %snps_hash;
  my $snps;

  if (!scalar(@snptypes)) {
    $snps = $slice->get_all_DnaAlignFeatures;
  } else {
    foreach my $snptype (@snptypes) {
      $snps_hash{$snptype} = $slice->get_all_DnaAlignFeatures($snptype);
      print "Got " . scalar(@{$snps_hash{$snptype}}) . " $snptype SNPs\n";
      push @$snps,@{$snps_hash{$snptype}};
    }
  }

  print "Done fetching SNPs\n";


  print "Starting sorting SNPs\n";
  my  @snps_array = sort {$a->start <=> $b->start} @$snps;
  $snps = \@snps_array;
  print "Done sorting SNPs\n";

# Go through the genes sorting the transcripts
  print "Starting sorting transcripts\n";
  foreach my $gene (@genes) {
    foreach my $trans (@{$gene->get_all_Transcripts}) {
      $trans->sort;
    }
  }
  print "Done sorting transcripts\n";

  my $count=0;
  my $ngene = scalar(@genes);
  my %snpcodingtype;
  my $codingtype = undef;

  print "Starting coding type calcs\n";

  GENESNP: foreach my $snp (@$snps) {

    $codingtype = undef;

    if (!($count%1000)) {
      print ".";
    }
    $count++;


    GENE: for (my $ind=0; $ind < $ngene; $ind++) {
      my $gene = $genes[$ind]; 
      if (!defined($gene)) {
        next GENE;
      }
      if ($snp->end >= $gene->start && 
          $snp->start <= $gene->end) {
        if (!defined($codingtype)) { $codingtype = 'intronic'; }

        foreach my $trans (@{$gene->get_all_Transcripts}) {
          if ($snp->end >= $trans->start && 
              $snp->start <= $trans->end) {
            foreach my $exon (@{$trans->get_all_Exons}) {
              if ($exon->overlaps($snp,'ignore')) {
                if ($trans->translation && 
                    $snp->end >= $trans->coding_region_start && 
                    $snp->start <= $trans->coding_region_end) {
                  $codingtype = 'coding';
           
                  next GENESNP;
                } else {
                  $codingtype = 'noncoding';
                }
              }
            }
          }
        }
      } elsif ($gene->start > $snp->end) {
        next GENESNP;
      } elsif ($gene->end < $snp->start) {
        #print "Removing gene " . $genes[$ind]->stable_id . "\n";
        $genes[$ind] = undef;
      }
    }
  } continue {
    if (defined($codingtype)) {
      $snpcodingtype{$snp->analysis->logic_name}{$codingtype}++;
    } else {
      $snpcodingtype{$snp->analysis->logic_name}{intergenic}++;
    }     
  }
  print "\n\nCODING TYPE COUNTS\n\n";
  foreach my $logic_name (keys %snpcodingtype) {
    print  " SNP Type $logic_name\n";
    foreach my $codingtype (keys %{$snpcodingtype{$logic_name}}) {
      print  "  Type $codingtype count = " . $snpcodingtype{$logic_name}{$codingtype} . "\n";
    }
    print  "\n";
  }
  print "Done coding type calcs\n";

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


