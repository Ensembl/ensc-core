#!/usr/bin/perl
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

#
# Script for adding new c source files to the jalview Makefile.
# SMJS 2/01
#
# NOTES: 
#   1. If you add new directories to the Jalview hierarchy this 
#      script will need changing.

use strict;

our ($opt_s,$opt_p);

use Getopt::Std;

use Env qw(ENSC_ROOT);

use File::PathConvert qw(rel2abs);

if (!defined ($ENSC_ROOT)) {
   die "ENSC_ROOT must be set\n";
}

my $SrcDir = "$ENSC_ROOT/src";

open MFILE, "<$SrcDir/Makefile" || die "Failed opening $SrcDir/Makefile\n";

my $DepLineStr = "#DO NOT REMOVE THIS LINE OR ADD ANYTHING AFTER IT - ADDMAKE DEPENDS ON IT FOR MAKING DEPENDENCIES";

system("cp $SrcDir/Makefile $SrcDir/Makefile.bck");
my $insect = 0;


getopts("ps:");

if ($#ARGV == -1) {
  usage();
}

my $section;
if ($opt_s ne "") {
  $section = $opt_s;
} else {
  $section = "libensc.a";
}

#read makefile to find section
my $line;
my @Before;
my @After;
my @Files;
while (defined($line = <MFILE>)) {
   chomp $line;
   if ($insect == 1) {
      if ( $line =~ /^$/) { 
         $insect = 2;
      } else {
         if ($line !~ /\\$/) {
           $insect = 4;
         }
         $line =~ s/\t//g;
         $line =~ s/ *\\ *$//;
         push @Files,$line;
      }
   } 
   if ($insect == 0) {
      push @Before,$line;
   }
   if ($line =~ /^$section\t:/) {
      $insect = 1;
   }
   if ($line =~ /$DepLineStr/) {
      $insect = 3;
   }
   if ($insect == 2) {
      push @After,$line;
   }
   if ($insect == 4) {
      $insect = 2;
   }
}

if ($insect != 3) {
   die "Failed parsing Makefile - it has not been modified.\n";
}

my @NewSrcFiles;
foreach my $newfile (@ARGV) {
   my $realname = "";
   #my $realname = `realpath $newfile`;
   my $cutname = $realname;

   if ($realname eq "") {
      $cutname = fakePath($newfile);
   }
   chomp $cutname;
   print $cutname . "\n";
   if (!($cutname =~ / *$SrcDir\//)) {
      die "New file $cutname didn't have $SrcDir in its path component\n";
   }
   if (!($cutname =~ /\.c$/)) {
      die "New file $cutname didn't have .c extension\n";
   }
   $cutname =~ s/ *$SrcDir\///;
   push @NewSrcFiles,$cutname;
   $cutname = toObjectFile($cutname);
   my $found = 0;
   foreach my $file (@Files) {
      if ($file eq $cutname) {
         $found = 1;
      }
   }
   if (!$found) {
      push @Files,$cutname;
   } else {
      print STDERR "$cutname already in Makefile\n";
      pop @NewSrcFiles;
   }
}

close MFILE;

print STDERR "Opening $SrcDir/Makefile for write\n";

open MFILE,">$SrcDir/Makefile";

foreach $line (@Before) {
   print  MFILE "$line\n";
}

if (! $opt_p) {
   my @pruned_files;
   foreach my $file (sort @Files) {
      my $srcfile = toAbsSrcFile($file,$SrcDir);
   
      if (! -e $srcfile ) {
         print STDERR "File $srcfile doesn't exist - pruning from Makefile\n";
      } else {
         push @pruned_files,$file;
      }
   }
   @Files = @pruned_files;
}


my $count = $#Files;
foreach my $file (sort @Files) {
   print MFILE "\t$file";
   my $srcfile = toAbsSrcFile($file,$SrcDir);
   if (! -e $srcfile ) {
      print STDERR "File $srcfile doesn't exist - added anyway\n";
   }
   if ($count) {
      print MFILE " \\";
   }  
   print MFILE "\n";
   $count--;
}

foreach $line (@After) {
   print MFILE "$line\n";
}

foreach my $file (@Files) {
   my $srcfile = $file;
   $srcfile =~ s/objectfiles\///;
   $srcfile =~ s/\.o$/.c/;
   push @NewSrcFiles,$srcfile;
}

print MFILE "$DepLineStr\n\n";

close MFILE;   
print STDERR "Finished writing to $SrcDir/Makefile\n";

sub toAbsSrcFile {
  my $objectname = shift;
  my $SrcDir = shift;

  return $SrcDir . "/" . toSrcFile($objectname);
}

sub toSrcFile {
  my $objectname = shift;
  my $srcname = $objectname;

  $srcname =~ s/^objectfiles\///;
  $srcname =~ s/\.o/.c/;
 
  return $srcname;
}
 
sub toObjectFile {
  my $srcname = shift;
  my $objectname = $srcname;

   $objectname =~ s/\.c/.o/;
   $objectname = "objectfiles/" . $objectname;

   return $objectname;
}

sub fakePath {
   my $file = shift;
   
#   return File::PathConvert::rel2abs($file)
   return rel2abs($file)
  
}

sub usage {
  print STDERR "Usage: addmake.pl [-s <section>] [-d] [-c] <files>\n";
  print STDERR "  Options:\n";
  print STDERR "        -s <section> - Specify which target to add new files to (defaults to gui)\n";
  die;
}
