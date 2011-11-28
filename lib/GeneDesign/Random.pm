package GeneDesign::Random;
require Exporter;

use GeneDesign::Codons(pattern_finder);

use strict;
use warnings;

use vars qw($VERSION @ISA @EXPORT_OK %EXPORT_TAGS);
$VERSION = 3.00;

@ISA = qw(Exporter);
@EXPORT_OK = qw(
  randDNA
);

%EXPORT_TAGS =  (all => [qw(randDNA)]);

################################################################################
###########################                          ###########################
################################################################################

sub randDNA
{
  my ($length, $ATperc, $stopswit, $CODON_TABLE) = @_;
  my $ATtotal = int( ( $ATperc * $length / 100 )  + .5) ;
  my $Acount  = int( rand( $ATtotal ) + .5 );
  my $Tcount  = $ATtotal - $Acount;
  my $GCtotal = $length - $ATtotal;
  my $Gcount  = int( rand( $GCtotal ) + .5 );
  my $Ccount  = $GCtotal - $Gcount;  
  my @randomarray = shuffle( split( '', ('A' x $Acount) . ('T' x $Tcount) . ('C' x $Ccount) . ('G' x $Gcount) ) );
  if ($stopswit == 2)
  {
    my $randDNA = join('', @randomarray);
    foreach (pattern_finder($randDNA, "*", 2, 1, $CODON_TABLE))
    {
      substr($randDNA, (($_+1)*3) - 3, 3) = scalar reverse(substr($randDNA, (($_+1)*3) - 3, 3));
      substr($randDNA, (($_+1)*3) - 2, 2) = scalar reverse(substr($randDNA, (($_+1)*3) - 2, 2)) if (int(rand(1)+.5) == 1);
    }
    return $randDNA;
  }
  else
  {
    return join("", @randomarray);
  }
}

1;

__END__

=head1 NAME

GeneDesign::Random

=head1 VERSION

Version 3.00

=head1 DESCRIPTION

  Random DNA Generators
  
=head1 Functions

=head2 randDNA()
  takes a target length and an AT percentage and generates a random nucleotide
  sequence, with or without stops in the first frame
  in: nucleotide sequence length (scalar), 
      AT percentage (0 <= scalar <= 100), 
      stop codon prevention(0 stops, 1 no stops), 
      codon table (hash reference)
  out: nucleotide sequence (string)
  
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
