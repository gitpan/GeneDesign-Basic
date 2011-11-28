package GeneDesign::SufTree;
use 5.006;
use strict;

## Modified from code by John Kloss

my %AA_KEYS =	('A'=> 0, 'C'=> 1, 'D'=> 2, 'E'=> 3, 'F'=> 4, 'G'=> 5, 'H'=> 6, 
               'I'=> 7, 'K'=> 8, 'L'=> 9, 'M'=>10, 'N'=>11, 'P'=>12, 'Q'=>13, 
               'R'=>14, 'S'=>15, 'T'=>16, 'V'=>17, 'W'=>18, 'Y'=>19, '*'=>20	);

################################################################################
###########################                          ###########################
################################################################################

sub new_aa
{
	my ($class) = @_;
	my $self = 
	{
		root => [ ( undef ) x 23 ]
	};
	bless $self, $class;
	return $self;
}

sub root 
{ 
	my ($self, $root) = @_;	
	$self->{ root} = $root if defined($root);	
	return $self->{root};
}

sub add_aa_paths
{
  my ($self, $ehshref) = @_;
  foreach my $peptide (keys %$ehshref)
  {
    my $next = $self->{ root };
    foreach my $aa (split '', $peptide)
    {
      my $charnum = $AA_KEYS{$aa};
      $next->[$charnum] ||= [ (undef)  x 23 ];
      $next = $next->[$charnum];
    }
    $next->[21] = $$ehshref{$peptide};
    $next->[22] = $peptide;
  }
}

sub find_aa_paths
{
  my ($self, $protein) = @_;
  my @locations;
  my @seq = split '', $protein;
  for my $seq_idx (0..scalar(@seq))
  {
    my $cur_idx = $seq_idx;
    my $ref_idx = $AA_KEYS{$seq[$seq_idx]};
    my $ref     = $self->{root};
    while (++$cur_idx < scalar(@seq) and $ref)
    {
      if ($ref->[21]) 
      {  
        push @locations, [$_, $seq_idx+1, $ref->[22]] foreach (@{$ref->[21]});
      }
      $ref_idx = $AA_KEYS{$seq[$cur_idx]};
      $ref = $ref->[$ref_idx];
    }
  }
  return @locations;
}

1;

__END__

=head1 NAME

GeneDesign::SufTree - A suffix tree implementation for restriction enzyme search

=head1 VERSION

Version 3.00

=head1 DESCRIPTION

  GeneDesign uses this object to parse peptide sequences for restriction enzyme
  recognition site possibilities

=head1 Functions

=head2 new_aa()

=head2 add_aa_paths()

=head2 find_aa_paths()

=head1 AUTHOR

Sarah Richardson <notadoctor@jhu.edu>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011, GeneDesign developers
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Johns Hopkins nor the
      names of the developers may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE DEVELOPERS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=cut
