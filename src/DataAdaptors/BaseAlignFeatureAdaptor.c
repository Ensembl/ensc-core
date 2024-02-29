/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "BaseAlignFeatureAdaptor.h"

sub fetch_all_by_RawContig_and_pid {
  my($self, $contig, $pid, $logic_name) = @_;

  my $constraint;

  #get the primary table alias
  my @tabs = $self->_tables;
  my $alias = $tabs[0]->[1];

  if(defined $pid) {
    $constraint = "${alias}.perc_ident > $pid";
  }

  return $self->fetch_all_by_RawContig_constraint($contig, 
						  $constraint, 
						  $logic_name);
}

sub fetch_all_by_Slice_and_pid {
  my ($self,$slice,$pid, $logic_name) = @_;
  my $constraint;


  #get the primary table alias
  my @tabs = $self->_tables;
  my $alias = $tabs[0]->[1];

  if(defined $pid) {
    $constraint = "${alias}.perc_ident > $pid";
  }

  if(defined $pid){
    $constraint = "perc_ident > $pid";
  }

  return $self->fetch_all_by_Slice_constraint($slice, $constraint, 
					      $logic_name);
}  
