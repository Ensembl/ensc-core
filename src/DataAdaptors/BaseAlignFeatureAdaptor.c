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
