package GeneDesign::Codons;
require Exporter;

use GeneDesign::Basic qw(regres compareseqs %AACIDS $ambnt %NTIDES complement);
use List::Util qw(shuffle first);

use strict;
use warnings;

use vars qw($VERSION @ISA @EXPORT_OK %EXPORT_TAGS);
$VERSION = 3.00;

@ISA = qw(Exporter);
@EXPORT_OK = qw(
  translate
  reverse_translate
  degcodon_to_aas
  amb_translation
  amb_transcription
  pattern_finder
  pattern_remover
  pattern_adder
  pattern_aligner
  change_codons
  rscu_parser
  define_codon_table
  define_reverse_codon_table
  define_RSCU_values 
  RSCU_filter
  define_codon_percentages
  index_codon_percentages
  codon_count
  generate_RSCU_values
  define_aa_defaults
  $strcodon
);

%EXPORT_TAGS =  (all => [qw(define_aa_defaults define_codon_table 
  define_reverse_codon_table define_RSCU_values RSCU_filter 
  define_codon_percentages index_codon_percentages pattern_remover pattern_adder
  pattern_aligner pattern_finder change_codons reverse_translate 
  amb_transcription amb_translation degcodon_to_aas translate codon_count 
  generate_RSCU_values rscu_parser $strcodon)]);

our $strcodon  = qr/[ATCG]{3}/;

################################################################################
######################## Codon and Amino Acid Functions ########################
################################################################################

sub translate
{
  my ($nucseq, $swit, $CODON_TABLE) = @_;
#  return if ( ! $nucseq  ||  $nucseq =~ $ambnt);
  $nucseq = complement($nucseq, 1) if ($swit < 0);
  my $peptide = "";
  for (my $offset = abs($swit)-1; $offset < length($nucseq); $offset += 3)
  {
    my $codon = substr($nucseq, $offset, 3);
    $peptide .= $$CODON_TABLE{$codon} if (exists $$CODON_TABLE{$codon});
  }
  return $peptide;
}

sub reverse_translate
{
  my($aaseq, $codonhash) = @_;
  my $newseq = "";
  $newseq .= $$codonhash{$_} foreach (split('', $aaseq));
  return $newseq;
}

sub degcodon_to_aas 
{
  my ($codon, $CODON_TABLE, $xlationref) = @_;  
  return if ( ! $codon  ||  length($codon) != 3);
  my (@answer, %temphash) = ((), ());
  if (exists $$xlationref{$codon})
  {
    return @{$$xlationref{$codon}};
  }
  elsif ($codon eq "NNN")
  {
    %temphash = map { $_ => 1} values %$CODON_TABLE;
    @answer = keys %temphash;
  }
  else
  {
    my $reg = regres($codon, 1);
    %temphash = map {$$CODON_TABLE{$_}  => 1} grep { $_ =~ $reg } keys %$CODON_TABLE;
    @answer = keys %temphash; 
  }
  $$xlationref{$codon} = \@answer;
  return @answer;
}

sub amb_translation
{
  my ($site, $CODON_TABLE, $swit, $xlationref) = @_;
  $site = 'NN' . $site if (!$swit);
  my (%RES, @SEED, @NEW);
  for (my $j = 0; $j < 3; $j++)
  {
    my $gothrough = 0;
    for (my $offset = $j; $offset < (int(length($site))); $offset +=3)
    {
      my $tempcodon = substr($site, $offset, 3);
      if (!$swit)
      {
        $tempcodon .= 'N' while (length($tempcodon) < 3);
      }
      if ($gothrough == 0)
      {
        @SEED = degcodon_to_aas($tempcodon, $CODON_TABLE, $xlationref) ;
      }
      else
      {
        @NEW  = degcodon_to_aas($tempcodon, $CODON_TABLE, $xlationref);
        @SEED = combine(\@SEED, \@NEW);
      }
      $gothrough++;
    }
    $RES{$_}++ foreach(@SEED);
  }
  return keys %RES;
}

sub amb_transcription
{
  my ($ntseq, $CODON_TABLE, $pepseq) = @_;
  my (@SEED, @NEW) = ((), ());
  my $offset = 0;
  if ( !$pepseq )
  {
    while ($offset < length($ntseq))
    {
      my $template = substr($ntseq, $offset, 3);
      my $regtemp = regres($template);
      if ($template !~ $ambnt)
      {
        @SEED = ( $template ) if ($offset == 0);
        @NEW  = ( $template ) if ($offset >  0);
      }
      else
      {
        my @TEMP = grep { $_ =~ $regtemp } keys %$CODON_TABLE;
        @SEED = @TEMP if ($offset == 0);
        @NEW  = @TEMP if ($offset >  0);
      }
      unless ($offset == 0) 
      {
        @SEED = combine(\@SEED, \@NEW);
      }
      $offset += 3;
    }
  }
  else
  {
    my $REV_CODON = define_reverse_codon_table($CODON_TABLE);
    my @each = split("", $pepseq);
    while ($offset < scalar(@each))
    {
      my $codon = substr($ntseq, $offset*3, 3);
      my $peptide = $each[$offset];
      my @TEMP = grep {$_ =~ regres($codon, 1)} @{$$REV_CODON{$peptide}};
      @SEED = @TEMP if ($offset == 0);
      @NEW  = @TEMP if ($offset >  0);
      unless ($offset == 0) 
      {
        @SEED = combine(\@SEED, \@NEW);
      }
      $offset++;
    }
  }
  my %hsh = map {$_ => 1 } @SEED;
  return keys %hsh;
#  return grep {translate($_, 1, $hashref) eq $pepseq} keys %SEED_TOTAL if ($pepseq);
}

#### combine ####
# meant to work with the amb_trans* functions. 
# Basically builds a list of tree nodes.
# in: 2 x peptide lists (array reference)
# out: combined list of peptides (vector)
sub combine
{
  my ($arr1ref, $arr2ref) = @_;
  my @arr3 = ();
  foreach my $do (@$arr1ref)
  {
    push @arr3, $do . $_ foreach (@$arr2ref)
  }
  return @arr3;
}

sub pattern_finder
{
  my ($strand, $pattern, $swit, $frame, $CODON_TABLE) = @_;
  my @positions = ();
  if ($swit == 2)
  {
    return if (! $frame || ! $CODON_TABLE);
    $strand = translate($strand, $frame, $CODON_TABLE)
  }
  my $exp = regres($pattern, $swit);
  while ($strand =~ /(?=$exp)/ig)
  {
    push @positions, (pos $strand);
  }
  return @positions;
}

sub pattern_remover
{
  my ($critseg, $pattern, $CODON_TABLE, $RSCU_TABLE) = @_;
  my @patterns = @{$pattern};
  my $REV_CODON_TABLE = define_reverse_codon_table($CODON_TABLE); 
  my %changes;
  # for each codon position, get RSCU difference possibilities
  for (my $offset = 0; $offset < length($critseg); $offset += 3)  
  {
    my $codono = substr($critseg, $offset, 3);
    foreach my $codonp (  grep { $codono ne $_} 
                @{$$REV_CODON_TABLE{$$CODON_TABLE{$codono}}})
    {
      my $id = $offset . "." . $codono . "." . $codonp;
      $changes{$id}->{DRSCU} = abs($$RSCU_TABLE{$codonp} - $$RSCU_TABLE{$codono});
      $changes{$id}->{OFFSET} = $offset;
      $changes{$id}->{OLDCOD} = $codono;
      $changes{$id}->{NEWCOD} = $codonp;
    }
  }
  #Replace the least different codons one at a time, take first solution
  foreach my $id (sort {$changes{$a}->{DRSCU} <=> $changes{$b}->{DRSCU}} keys %changes)
  {
    my $copy = $critseg;
    substr($copy, $changes{$id}->{OFFSET}, 3) = $changes{$id}->{NEWCOD};
    next if ($copy =~ $patterns[0] || complement($copy, 1) =~ $patterns[0]);
    next if ($patterns[1] && ($copy =~ $patterns[1] || complement($copy, 1) =~ $patterns[1]));
    return $copy;
  }
  #Try pairwise combinations of codons, sorted by sum of RSCU difference
  my @singles = keys %changes;
  my @pairs;
  for my $x (0..scalar(@singles)-1)
  {
    for my $y ($x+1..scalar(@singles)-1)
    {
      if ($changes{$singles[$x]}->{OFFSET} != $changes{$singles[$y]}->{OFFSET})
      {
        push @pairs, [$singles[$x], $singles[$y]];
      }
    }
  }
  @pairs = sort {($changes{$a->[0]}->{DRSCU} + $changes{$a->[1]}->{DRSCU})
            <=>  ($changes{$b->[0]}->{DRSCU} + $changes{$b->[1]}->{DRSCU})} @pairs;
  foreach my $pair (@pairs)
  {
    my ($one, $two) = @$pair;
    my $copy = $critseg;
    substr($copy, $changes{$one}->{OFFSET}, 3) = $changes{$one}->{NEWCOD};
    substr($copy, $changes{$two}->{OFFSET}, 3) = $changes{$two}->{NEWCOD};
    next if ($copy =~ $patterns[0] || complement($copy, 1) =~ $patterns[0]);
    next if ($patterns[1] && ($copy =~ $patterns[1] || complement($copy, 1) =~ $patterns[1]));
    return $copy;
  }  
  return 0;#$critseg;
}

sub pattern_adder  
{
  my ($oldpatt, $newpatt, $CODON_TABLE, $xlationref) = @_;
  #assume that critseg and pattern come in as complete codons 
  # (i.e., have been run through pattern_aligner)
  my $REV_CODON_TABLE = define_reverse_codon_table($CODON_TABLE);
  my $copy = "";
  for (my $offset = 0; $offset < length($oldpatt); $offset += 3)
  {
    my $curcod = substr($oldpatt, $offset, 3);
    my $curtar = substr($newpatt, $offset, 3);
    foreach my $g (degcodon_to_aas($curcod, $CODON_TABLE, $xlationref))
    {
      $copy .= $curcod =~ regres($curtar)  
           ? $curcod
           : first { compareseqs($curtar, $_) } @{$$REV_CODON_TABLE{$g}};
      # print "\t\tpatadd\t($curcod, $curtar)\t$copy<br>\n";
    }
  }
  # print "\t\tpatadd $copy from $oldpatt\n";
  return length($copy) == length($oldpatt)  ?  $copy  :  0;
}

sub pattern_aligner
{
  my ($critseg, $pattern, $peptide, $CODON_TABLE, $swit, $xlationref) = @_;
  $swit = 0 if (!$swit);
  my $diff = length($critseg) - length($pattern);
  my ($newpatt, $nstring, $rounds, $offset, $check, $pcheck) = ("", "N" x $diff, 0, 0, "", "");
  #  print "seeking $pattern for $peptide from $critseg...\n";
  while ($check ne $peptide && $rounds <= $diff*2 + 1)
  {
    $newpatt = $rounds <= $diff  
      ?  substr($nstring, 0, $rounds) . $pattern  
      :  substr($nstring, 0, ($rounds-3)) . complement($pattern, 1);
    $newpatt .=  "N" while (length($newpatt) != length($critseg));
    #  print "\t$newpatt\n";
    my ($noff, $poff) = (0, 0);
    $check = "";
    while ($poff < length($peptide))
    {
      my @possibles = degcodon_to_aas( substr($newpatt, $noff, 3), $CODON_TABLE, $xlationref );
      #   print "\t\t@possibles\n";
      $check .= $_ foreach( grep { substr($peptide, $poff, 1) eq $_ } @possibles);
      $noff += 3;
      $poff ++;
    }
    $pcheck = translate(substr($critseg, $offset, length($peptide) * 3), 1, $CODON_TABLE);
    #      print "\t\t$check, $pcheck, $offset\n";
    $rounds++;
    $offset += 3 if ($rounds % 3 == 0);
#    $check = "" if ( $pcheck !~ $check);
  }
  $newpatt = "0" if ($check ne $peptide);
  # print "\t\tpataln $check, $pcheck, $rounds, $newpatt\n" if ($check ne $peptide);
  return $swit == 1  ?  ($newpatt, $rounds-1)  :  $newpatt;
}

sub change_codons
{
  my ($oldseq, $CODON_TABLE, $RSCU_VALUES, $swit, $tag) = @_;
  my $REV_CODON_TABLE = define_reverse_codon_table($CODON_TABLE);  
  my ($offset, $newcod, $curcod, $newseq, $aa) = (0, undef, undef, undef, undef);
  $tag = 0 if (! $tag);
  while ($offset < length($oldseq))
  {
    $curcod = substr($oldseq, $offset, 3);
    $newcod = $curcod;
    $aa = $$CODON_TABLE{$curcod};
    my @posarr = sort { $$RSCU_VALUES{$b} <=> $$RSCU_VALUES{$a} } 
           grep { exists($$RSCU_VALUES{$_}) }
           @{$$REV_CODON_TABLE{$aa}};
    if (scalar(@posarr) != 1 && $aa ne '*')
    {  
      if    ($swit == 0)  #Random
      {
        $newcod = first { $_ ne $curcod } shuffle @posarr;
      }
      elsif ($swit == 1)  #Optimal
      {
        $newcod = first {    1    } @posarr;
      }
      elsif ($swit == 2)  #Less Optimal
      {
        $newcod = first { $_ ne $curcod } @posarr;
      }
      elsif ($swit == 3)  #Most Different
      {
        my $lastbase = substr $curcod, 2, 1;  
        my $frstbase = substr $curcod, 0, 1;
        if ($tag && $offset == 0)
        {
          $newcod = first {  substr($_, 0, 2) eq substr($curcod, 0, 2) 
                  &&  substr($_, 2, 1) =~ sitver($lastbase, 1)}
                @posarr;
          if (!$newcod)
          {
            $newcod = first {  substr($_, 0, 2) eq substr($curcod, 0, 2) 
                    &&  substr($_, 2, 1) ne $lastbase}
                  @posarr;
          }
          if (!$newcod)
          {
            die ("can't find a better cod for $curcod $aa<br>");
          }
        }
        elsif  (scalar(@posarr) == 2)    
        {  
          $newcod = first { $_ ne $curcod } @posarr;  
        }
        elsif  (scalar(@posarr) <= 4)
        {
          $newcod = first { (substr $_, 2, 1) =~ sitver($lastbase, 1) } @posarr;  
          if (!$newcod)
          {
            $newcod = first { (substr $_, 2) ne $lastbase } @posarr;
          }
        }
        elsif  (scalar(@posarr) > 4)
        {  
          $newcod = first {  ((substr $_, 2) !~ sitver($lastbase, 0) ) 
                  && ((substr $_, 0, 1) ne $frstbase)} 
                @posarr;
          if (!$newcod)
          {
            $newcod = first {  ((substr $_, 2) ne $lastbase ) 
                    && ((substr $_, 0, 1) ne $frstbase)} 
                  @posarr;
          }
        }
      #  print "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$tag, $curcod, $newcod !<br><Br>";
      }
      elsif ($swit == 4)  #Least Different
      {
        my @sorarr = sort {abs($$RSCU_VALUES{$a} - $$RSCU_VALUES{$curcod}) <=> abs($$RSCU_VALUES{$b} - $$RSCU_VALUES{$curcod})} @posarr;
        $newcod = first { $_ ne $curcod } @sorarr;
        $newcod = $curcod if (abs($$RSCU_VALUES{$newcod} - $$RSCU_VALUES{$curcod}) > 1);
      }
    }
    $newseq .= $newcod;
    $offset += 3;
  }
  return $newseq;
}

sub sitver
{
  my ($base, $swit) = @_;
  $swit = 0 if (!$swit);
  if ($swit == 1) #return set of transversions
  {
    return $base =~ $NTIDES{Y}  ?  $NTIDES{R}  :  $NTIDES{Y};
  }
  else      #return set of transitions
  {
    return $base =~ $NTIDES{Y}  ?  $NTIDES{Y}  :  $NTIDES{R};
  }
}

sub rscu_parser
{
  my ($instr) = @_;
  my $rscuhsh = {};
  foreach my $pre (split(/[\n\r]/, $instr))
  {
    my @trip = split(/ /, $pre);
    my $id = $trip[0];
    my $val = $trip[2];
    $$rscuhsh{$id} = $val;
  }
  return $rscuhsh;
}

sub define_codon_table
{
  my ($swit) = @_;
  return if (!$swit);
  my $CODON_TABLE = {};
  if ($swit != 8)
  {
    $$CODON_TABLE{"TTT"} = "F";  $$CODON_TABLE{"TTC"} = "F";  $$CODON_TABLE{"TTA"} = "L";  $$CODON_TABLE{"TTG"} = "L";  
    $$CODON_TABLE{"CTT"} = "L";  $$CODON_TABLE{"CTC"} = "L";  $$CODON_TABLE{"CTA"} = "L";  $$CODON_TABLE{"CTG"} = "L";
    $$CODON_TABLE{"ATT"} = "I";  $$CODON_TABLE{"ATC"} = "I";  $$CODON_TABLE{"ATA"} = "I";  $$CODON_TABLE{"ATG"} = "M";  
    $$CODON_TABLE{"GTT"} = "V";  $$CODON_TABLE{"GTC"} = "V";  $$CODON_TABLE{"GTA"} = "V";  $$CODON_TABLE{"GTG"} = "V";  
    $$CODON_TABLE{"TCT"} = "S";  $$CODON_TABLE{"TCC"} = "S";  $$CODON_TABLE{"TCA"} = "S";  $$CODON_TABLE{"TCG"} = "S";  
    $$CODON_TABLE{"CCT"} = "P";  $$CODON_TABLE{"CCC"} = "P";  $$CODON_TABLE{"CCA"} = "P";  $$CODON_TABLE{"CCG"} = "P";
    $$CODON_TABLE{"ACT"} = "T";  $$CODON_TABLE{"ACC"} = "T";  $$CODON_TABLE{"ACA"} = "T";  $$CODON_TABLE{"ACG"} = "T";  
    $$CODON_TABLE{"GCT"} = "A";  $$CODON_TABLE{"GCC"} = "A";  $$CODON_TABLE{"GCA"} = "A";  $$CODON_TABLE{"GCG"} = "A";  
    $$CODON_TABLE{"TAT"} = "Y";  $$CODON_TABLE{"TAC"} = "Y";  $$CODON_TABLE{"TAA"} = "*";  $$CODON_TABLE{"TAG"} = "*";  
    $$CODON_TABLE{"CAT"} = "H";  $$CODON_TABLE{"CAC"} = "H";  $$CODON_TABLE{"CAA"} = "Q";  $$CODON_TABLE{"CAG"} = "Q";
    $$CODON_TABLE{"AAT"} = "N";  $$CODON_TABLE{"AAC"} = "N";  $$CODON_TABLE{"AAA"} = "K";  $$CODON_TABLE{"AAG"} = "K";  
    $$CODON_TABLE{"GAT"} = "D";  $$CODON_TABLE{"GAC"} = "D";  $$CODON_TABLE{"GAA"} = "E";  $$CODON_TABLE{"GAG"} = "E";  
    $$CODON_TABLE{"TGT"} = "C";  $$CODON_TABLE{"TGC"} = "C";  $$CODON_TABLE{"TGA"} = "*";  $$CODON_TABLE{"TGG"} = "W";  
    $$CODON_TABLE{"CGT"} = "R";  $$CODON_TABLE{"CGC"} = "R";  $$CODON_TABLE{"CGA"} = "R";  $$CODON_TABLE{"CGG"} = "R";  
    $$CODON_TABLE{"AGT"} = "S";  $$CODON_TABLE{"AGC"} = "S";  $$CODON_TABLE{"AGA"} = "R";  $$CODON_TABLE{"AGG"} = "R";  
    $$CODON_TABLE{"GGT"} = "G";  $$CODON_TABLE{"GGC"} = "G";  $$CODON_TABLE{"GGA"} = "G";  $$CODON_TABLE{"GGG"} = "G";
  }
  elsif ($swit == 8)
  {
    $$CODON_TABLE{"TTT"} = "F";  $$CODON_TABLE{"TTC"} = "F";  $$CODON_TABLE{"TTA"} = "L";  $$CODON_TABLE{"TTG"} = "L";  
    $$CODON_TABLE{"CTT"} = "L";  $$CODON_TABLE{"CTC"} = "L";  $$CODON_TABLE{"CTA"} = "L";  $$CODON_TABLE{"CTG"} = "L";
    $$CODON_TABLE{"ATT"} = "I";  $$CODON_TABLE{"ATC"} = "I";  $$CODON_TABLE{"ATA"} = "I";  $$CODON_TABLE{"ATG"} = "M";  
    $$CODON_TABLE{"GTT"} = "V";  $$CODON_TABLE{"GTC"} = "V";  $$CODON_TABLE{"GTA"} = "V";  $$CODON_TABLE{"GTG"} = "V";  
    $$CODON_TABLE{"TCT"} = "S";  $$CODON_TABLE{"TCC"} = "S";  $$CODON_TABLE{"TCA"} = "S";  $$CODON_TABLE{"TCG"} = "S";  
    $$CODON_TABLE{"CCT"} = "P";  $$CODON_TABLE{"CCC"} = "P";  $$CODON_TABLE{"CCA"} = "P";  $$CODON_TABLE{"CCG"} = "P";
    $$CODON_TABLE{"ACT"} = "T";  $$CODON_TABLE{"ACC"} = "T";  $$CODON_TABLE{"ACA"} = "T";  $$CODON_TABLE{"ACG"} = "T";  
    $$CODON_TABLE{"GCT"} = "A";  $$CODON_TABLE{"GCC"} = "A";  $$CODON_TABLE{"GCA"} = "A";  $$CODON_TABLE{"GCG"} = "A";  
    $$CODON_TABLE{"TAT"} = "Y";  $$CODON_TABLE{"TAC"} = "Y";  $$CODON_TABLE{"TAA"} = "*";  $$CODON_TABLE{"TAG"} = "*";  
    $$CODON_TABLE{"CAT"} = "H";  $$CODON_TABLE{"CAC"} = "H";  $$CODON_TABLE{"CAA"} = "Q";  $$CODON_TABLE{"CAG"} = "Q";
    $$CODON_TABLE{"AAT"} = "N";  $$CODON_TABLE{"AAC"} = "N";  $$CODON_TABLE{"AAA"} = "K";  $$CODON_TABLE{"AAG"} = "K";  
    $$CODON_TABLE{"GAT"} = "D";  $$CODON_TABLE{"GAC"} = "D";  $$CODON_TABLE{"GAA"} = "E";  $$CODON_TABLE{"GAG"} = "E";  
    $$CODON_TABLE{"TGT"} = "C";  $$CODON_TABLE{"TGC"} = "C";  $$CODON_TABLE{"TGA"} = "W";  $$CODON_TABLE{"TGG"} = "W";  
    $$CODON_TABLE{"CGT"} = "R";  $$CODON_TABLE{"CGC"} = "R";  $$CODON_TABLE{"CGA"} = "R";  $$CODON_TABLE{"CGG"} = "R";  
    $$CODON_TABLE{"AGT"} = "S";  $$CODON_TABLE{"AGC"} = "S";  $$CODON_TABLE{"AGA"} = "R";  $$CODON_TABLE{"AGG"} = "R";  
    $$CODON_TABLE{"GGT"} = "G";  $$CODON_TABLE{"GGC"} = "G";  $$CODON_TABLE{"GGA"} = "G";  $$CODON_TABLE{"GGG"} = "G";
  }
  return $CODON_TABLE;
}

sub define_reverse_codon_table
{
  my ($CODON_TABLE) = @_;
  my $REV_CODON_TABLE = {};
  foreach my $codon (keys %$CODON_TABLE)
  {
    my $aa = $$CODON_TABLE{$codon};
    $$REV_CODON_TABLE{$aa} = [] if ( ! exists $$REV_CODON_TABLE{$aa} );
    push @{$$REV_CODON_TABLE{$aa}}, $codon;
  }
  return $REV_CODON_TABLE;
}

sub define_RSCU_values
{
  my ($swit) = @_;
  return if (!$swit);
  my %RSCU_TABLE;
  if ($swit == 1)    #Saccharomyces cerevisiae
  {
    $RSCU_TABLE{"TTT"} = 0.19;  $RSCU_TABLE{"TTC"} = 1.81;  $RSCU_TABLE{"TTA"} = 0.49;  $RSCU_TABLE{"TTG"} = 5.34;  
    $RSCU_TABLE{"CTT"} = 0.02;  $RSCU_TABLE{"CTC"} = 0.00;  $RSCU_TABLE{"CTA"} = 0.15;  $RSCU_TABLE{"CTG"} = 0.02;
    $RSCU_TABLE{"ATT"} = 1.26;  $RSCU_TABLE{"ATC"} = 1.74;  $RSCU_TABLE{"ATA"} = 0.00;  $RSCU_TABLE{"ATG"} = 1.00;  
    $RSCU_TABLE{"GTT"} = 2.07;  $RSCU_TABLE{"GTC"} = 1.91;  $RSCU_TABLE{"GTA"} = 0.00;  $RSCU_TABLE{"GTG"} = 0.02;  
    $RSCU_TABLE{"TCT"} = 3.26;  $RSCU_TABLE{"TCC"} = 2.42;  $RSCU_TABLE{"TCA"} = 0.08;  $RSCU_TABLE{"TCG"} = 0.02;  
    $RSCU_TABLE{"CCT"} = 0.21;  $RSCU_TABLE{"CCC"} = 0.02;  $RSCU_TABLE{"CCA"} = 3.77;  $RSCU_TABLE{"CCG"} = 0.00;
    $RSCU_TABLE{"ACT"} = 1.83;  $RSCU_TABLE{"ACC"} = 2.15;  $RSCU_TABLE{"ACA"} = 0.00;  $RSCU_TABLE{"ACG"} = 0.01;  
    $RSCU_TABLE{"GCT"} = 3.09;  $RSCU_TABLE{"GCC"} = 0.89;  $RSCU_TABLE{"GCA"} = 0.03;  $RSCU_TABLE{"GCG"} = 0.00;  
    $RSCU_TABLE{"TAT"} = 0.06;  $RSCU_TABLE{"TAC"} = 1.94;  $RSCU_TABLE{"TAA"} = 1.00;  $RSCU_TABLE{"TAG"} = 0.00;  
    $RSCU_TABLE{"CAT"} = 0.32;  $RSCU_TABLE{"CAC"} = 1.68;  $RSCU_TABLE{"CAA"} = 1.98;  $RSCU_TABLE{"CAG"} = 0.02;
    $RSCU_TABLE{"AAT"} = 0.06;  $RSCU_TABLE{"AAC"} = 1.94;  $RSCU_TABLE{"AAA"} = 0.16;  $RSCU_TABLE{"AAG"} = 1.84;  
    $RSCU_TABLE{"GAT"} = 0.70;  $RSCU_TABLE{"GAC"} = 1.30;  $RSCU_TABLE{"GAA"} = 1.98;  $RSCU_TABLE{"GAG"} = 0.02;  
    $RSCU_TABLE{"TGT"} = 1.80;  $RSCU_TABLE{"TGC"} = 0.20;  $RSCU_TABLE{"TGA"} = 0.00;  $RSCU_TABLE{"TGG"} = 1.00;  
    $RSCU_TABLE{"CGT"} = 0.63;  $RSCU_TABLE{"CGC"} = 0.00;  $RSCU_TABLE{"CGA"} = 0.00;  $RSCU_TABLE{"CGG"} = 0.00;  
    $RSCU_TABLE{"AGT"} = 0.06;  $RSCU_TABLE{"AGC"} = 0.16;  $RSCU_TABLE{"AGA"} = 5.37;  $RSCU_TABLE{"AGG"} = 0.00;  
    $RSCU_TABLE{"GGT"} = 3.92;  $RSCU_TABLE{"GGC"} = 0.06;  $RSCU_TABLE{"GGA"} = 0.00;  $RSCU_TABLE{"GGG"} = 0.02;
  }
  elsif ($swit == 2)  #Escherichia coli
  {
    $RSCU_TABLE{"TTT"} = 0.34;  $RSCU_TABLE{"TTC"} = 1.66;  $RSCU_TABLE{"TTA"} = 0.06;  $RSCU_TABLE{"TTG"} = 0.07;  
    $RSCU_TABLE{"CTT"} = 0.13;  $RSCU_TABLE{"CTC"} = 0.17;  $RSCU_TABLE{"CTA"} = 0.04;  $RSCU_TABLE{"CTG"} = 5.54;
    $RSCU_TABLE{"ATT"} = 0.48;  $RSCU_TABLE{"ATC"} = 2.51;  $RSCU_TABLE{"ATA"} = 0.01;  $RSCU_TABLE{"ATG"} = 1.00;  
    $RSCU_TABLE{"GTT"} = 2.41;  $RSCU_TABLE{"GTC"} = 0.08;  $RSCU_TABLE{"GTA"} = 1.12;  $RSCU_TABLE{"GTG"} = 0.40;  
    $RSCU_TABLE{"TCT"} = 2.81;  $RSCU_TABLE{"TCC"} = 2.07;  $RSCU_TABLE{"TCA"} = 0.06;  $RSCU_TABLE{"TCG"} = 0.00;  
    $RSCU_TABLE{"CCT"} = 0.15;  $RSCU_TABLE{"CCC"} = 0.02;  $RSCU_TABLE{"CCA"} = 0.42;  $RSCU_TABLE{"CCG"} = 3.41;
    $RSCU_TABLE{"ACT"} = 1.87;  $RSCU_TABLE{"ACC"} = 1.91;  $RSCU_TABLE{"ACA"} = 0.10;  $RSCU_TABLE{"ACG"} = 0.12;  
    $RSCU_TABLE{"GCT"} = 2.02;  $RSCU_TABLE{"GCC"} = 0.18;  $RSCU_TABLE{"GCA"} = 1.09;  $RSCU_TABLE{"GCG"} = 0.71;  
    $RSCU_TABLE{"TAT"} = 0.38;  $RSCU_TABLE{"TAC"} = 1.63;  $RSCU_TABLE{"TAA"} = 0.00;  $RSCU_TABLE{"TAG"} = 0.00;  
    $RSCU_TABLE{"CAT"} = 0.45;  $RSCU_TABLE{"CAC"} = 1.55;  $RSCU_TABLE{"CAA"} = 0.12;  $RSCU_TABLE{"CAG"} = 1.88;
    $RSCU_TABLE{"AAT"} = 0.02;  $RSCU_TABLE{"AAC"} = 1.98;  $RSCU_TABLE{"AAA"} = 1.63;  $RSCU_TABLE{"AAG"} = 0.37;  
    $RSCU_TABLE{"GAT"} = 0.51;  $RSCU_TABLE{"GAC"} = 1.49;  $RSCU_TABLE{"GAA"} = 1.64;  $RSCU_TABLE{"GAG"} = 0.36;  
    $RSCU_TABLE{"TGT"} = 0.60;  $RSCU_TABLE{"TGC"} = 1.40;  $RSCU_TABLE{"TGA"} = 0.00;  $RSCU_TABLE{"TGG"} = 1.00;  
    $RSCU_TABLE{"CGT"} = 4.47;  $RSCU_TABLE{"CGC"} = 1.53;  $RSCU_TABLE{"CGA"} = 0.00;  $RSCU_TABLE{"CGG"} = 0.00;  
    $RSCU_TABLE{"AGT"} = 0.13;  $RSCU_TABLE{"AGC"} = 0.93;  $RSCU_TABLE{"AGA"} = 0.00;  $RSCU_TABLE{"AGG"} = 0.00;  
    $RSCU_TABLE{"GGT"} = 2.27;  $RSCU_TABLE{"GGC"} = 1.68;  $RSCU_TABLE{"GGA"} = 0.00;  $RSCU_TABLE{"GGG"} = 0.04;
  }
  elsif ($swit == 3)  #H.sapiens
  {
    $RSCU_TABLE{"TTT"} = 0.27;  $RSCU_TABLE{"TTC"} = 1.73;  $RSCU_TABLE{"TTA"} = 0.05;  $RSCU_TABLE{"TTG"} = 0.31;  
    $RSCU_TABLE{"CTT"} = 0.20;  $RSCU_TABLE{"CTC"} = 1.42;  $RSCU_TABLE{"CTA"} = 0.15;  $RSCU_TABLE{"CTG"} = 3.88;
    $RSCU_TABLE{"ATT"} = 0.45;  $RSCU_TABLE{"ATC"} = 2.43;  $RSCU_TABLE{"ATA"} = 0.12;  $RSCU_TABLE{"ATG"} = 1.00;  
    $RSCU_TABLE{"GTT"} = 0.09;  $RSCU_TABLE{"GTC"} = 1.03;  $RSCU_TABLE{"GTA"} = 0.11;  $RSCU_TABLE{"GTG"} = 2.78;  
    $RSCU_TABLE{"TCT"} = 0.45;  $RSCU_TABLE{"TCC"} = 2.09;  $RSCU_TABLE{"TCA"} = 0.26;  $RSCU_TABLE{"TCG"} = 0.68;  
    $RSCU_TABLE{"CCT"} = 0.58;  $RSCU_TABLE{"CCC"} = 2.02;  $RSCU_TABLE{"CCA"} = 0.36;  $RSCU_TABLE{"CCG"} = 1.04;
    $RSCU_TABLE{"ACT"} = 0.36;  $RSCU_TABLE{"ACC"} = 2.37;  $RSCU_TABLE{"ACA"} = 0.36;  $RSCU_TABLE{"ACG"} = 0.92;  
    $RSCU_TABLE{"GCT"} = 0.45;  $RSCU_TABLE{"GCC"} = 2.38;  $RSCU_TABLE{"GCA"} = 0.36;  $RSCU_TABLE{"GCG"} = 0.82;
    $RSCU_TABLE{"TAT"} = 0.34;  $RSCU_TABLE{"TAC"} = 1.66;  $RSCU_TABLE{"TAA"} = 0.00;  $RSCU_TABLE{"TAG"} = 0.00;  
    $RSCU_TABLE{"CAT"} = 0.30;  $RSCU_TABLE{"CAC"} = 1.70;  $RSCU_TABLE{"CAA"} = 0.21;  $RSCU_TABLE{"CAG"} = 1.79;
    $RSCU_TABLE{"AAT"} = 0.33;  $RSCU_TABLE{"AAC"} = 1.67;  $RSCU_TABLE{"AAA"} = 0.34;  $RSCU_TABLE{"AAG"} = 1.66;  
    $RSCU_TABLE{"GAT"} = 0.36;  $RSCU_TABLE{"GAC"} = 1.64;  $RSCU_TABLE{"GAA"} = 0.26;  $RSCU_TABLE{"GAG"} = 1.74;
    $RSCU_TABLE{"TGT"} = 0.42;  $RSCU_TABLE{"TGC"} = 1.58;  $RSCU_TABLE{"TGA"} = 0.00;  $RSCU_TABLE{"TGG"} = 1.00;
    $RSCU_TABLE{"CGT"} = 0.38;  $RSCU_TABLE{"CGC"} = 2.72;  $RSCU_TABLE{"CGA"} = 0.31;  $RSCU_TABLE{"CGG"} = 1.53;  
    $RSCU_TABLE{"AGT"} = 0.31;  $RSCU_TABLE{"AGC"} = 2.22;  $RSCU_TABLE{"AGA"} = 0.22;  $RSCU_TABLE{"AGG"} = 0.84;  
    $RSCU_TABLE{"GGT"} = 0.34;  $RSCU_TABLE{"GGC"} = 2.32;  $RSCU_TABLE{"GGA"} = 0.29;  $RSCU_TABLE{"GGG"} = 1.05;
  }
  elsif ($swit == 4)  #Caenorhabditis elegans
  {
    $RSCU_TABLE{"TTT"} = 0.72;  $RSCU_TABLE{"TTC"} = 1.28;  $RSCU_TABLE{"TTA"} = 0.45;  $RSCU_TABLE{"TTG"} = 1.35;  
    $RSCU_TABLE{"CTT"} = 1.86;  $RSCU_TABLE{"CTC"} = 1.38;  $RSCU_TABLE{"CTA"} = 0.34;  $RSCU_TABLE{"CTG"} = 0.63;
    $RSCU_TABLE{"ATT"} = 1.52;  $RSCU_TABLE{"ATC"} = 1.23;  $RSCU_TABLE{"ATA"} = 0.25;  $RSCU_TABLE{"ATG"} = 1.00;  
    $RSCU_TABLE{"GTT"} = 1.67;  $RSCU_TABLE{"GTC"} = 1.11;  $RSCU_TABLE{"GTA"} = 0.52;  $RSCU_TABLE{"GTG"} = 0.70;
    $RSCU_TABLE{"TCT"} = 1.47;  $RSCU_TABLE{"TCC"} = 0.98;  $RSCU_TABLE{"TCA"} = 1.44;  $RSCU_TABLE{"TCG"} = 0.83;  
    $RSCU_TABLE{"CCT"} = 0.52;  $RSCU_TABLE{"CCC"} = 0.23;  $RSCU_TABLE{"CCA"} = 2.75;  $RSCU_TABLE{"CCG"} = 0.51;
    $RSCU_TABLE{"ACT"} = 1.34;  $RSCU_TABLE{"ACC"} = 1.02;  $RSCU_TABLE{"ACA"} = 1.15;  $RSCU_TABLE{"ACG"} = 0.49;  
    $RSCU_TABLE{"GCT"} = 1.64;  $RSCU_TABLE{"GCC"} = 1.06;  $RSCU_TABLE{"GCA"} = 0.99;  $RSCU_TABLE{"GCG"} = 0.31;
    $RSCU_TABLE{"TAT"} = 0.97;  $RSCU_TABLE{"TAC"} = 1.03;  $RSCU_TABLE{"TAA"} = 0.00;  $RSCU_TABLE{"TAG"} = 0.00;  
    $RSCU_TABLE{"CAT"} = 1.13;  $RSCU_TABLE{"CAC"} = 0.87;  $RSCU_TABLE{"CAA"} = 1.39;  $RSCU_TABLE{"CAG"} = 0.61;
    $RSCU_TABLE{"AAT"} = 1.10;  $RSCU_TABLE{"AAC"} = 0.90;  $RSCU_TABLE{"AAA"} = 0.84;  $RSCU_TABLE{"AAG"} = 1.16;  
    $RSCU_TABLE{"GAT"} = 1.36;  $RSCU_TABLE{"GAC"} = 0.64;  $RSCU_TABLE{"GAA"} = 1.15;  $RSCU_TABLE{"GAG"} = 0.85;      
    $RSCU_TABLE{"TGT"} = 1.14;  $RSCU_TABLE{"TGC"} = 0.86;  $RSCU_TABLE{"TGA"} = 0.00;  $RSCU_TABLE{"TGG"} = 1.00;  
    $RSCU_TABLE{"CGT"} = 1.84;  $RSCU_TABLE{"CGC"} = 0.73;  $RSCU_TABLE{"CGA"} = 1.07;  $RSCU_TABLE{"CGG"} = 0.31;  
    $RSCU_TABLE{"AGT"} = 0.76;  $RSCU_TABLE{"AGC"} = 0.52;  $RSCU_TABLE{"AGA"} = 1.79;  $RSCU_TABLE{"AGG"} = 0.26;  
    $RSCU_TABLE{"GGT"} = 0.70;  $RSCU_TABLE{"GGC"} = 0.28;  $RSCU_TABLE{"GGA"} = 2.85;  $RSCU_TABLE{"GGG"} = 0.16;
  }
  elsif ($swit == 5)  #Drosophila melanogaster
  {
    $RSCU_TABLE{"TTT"} = 0.12;  $RSCU_TABLE{"TTC"} = 1.88;  $RSCU_TABLE{"TTA"} = 0.03;  $RSCU_TABLE{"TTG"} = 0.69;  
    $RSCU_TABLE{"CTT"} = 0.25;  $RSCU_TABLE{"CTC"} = 0.72;  $RSCU_TABLE{"CTA"} = 0.06;  $RSCU_TABLE{"CTG"} = 4.25;
    $RSCU_TABLE{"ATT"} = 0.74;  $RSCU_TABLE{"ATC"} = 2.26;  $RSCU_TABLE{"ATA"} = 0.00;  $RSCU_TABLE{"ATG"} = 1.00;  
    $RSCU_TABLE{"GTT"} = 0.56;  $RSCU_TABLE{"GTC"} = 1.59;  $RSCU_TABLE{"GTA"} = 0.06;  $RSCU_TABLE{"GTG"} = 1.79;
    $RSCU_TABLE{"TCT"} = 0.87;  $RSCU_TABLE{"TCC"} = 2.74;  $RSCU_TABLE{"TCA"} = 0.04;  $RSCU_TABLE{"TCG"} = 1.17;  
    $RSCU_TABLE{"CCT"} = 0.42;  $RSCU_TABLE{"CCC"} = 2.73;  $RSCU_TABLE{"CCA"} = 0.62;  $RSCU_TABLE{"CCG"} = 0.23;
    $RSCU_TABLE{"ACT"} = 0.65;  $RSCU_TABLE{"ACC"} = 3.04;  $RSCU_TABLE{"ACA"} = 0.10;  $RSCU_TABLE{"ACG"} = 0.21;  
    $RSCU_TABLE{"GCT"} = 0.95;  $RSCU_TABLE{"GCC"} = 2.82;  $RSCU_TABLE{"GCA"} = 0.09;  $RSCU_TABLE{"GCG"} = 0.14;
    $RSCU_TABLE{"TAT"} = 0.23;  $RSCU_TABLE{"TAC"} = 1.77;  $RSCU_TABLE{"TAA"} = 0.00;  $RSCU_TABLE{"TAG"} = 0.00;  
    $RSCU_TABLE{"CAT"} = 0.29;  $RSCU_TABLE{"CAC"} = 1.71;  $RSCU_TABLE{"CAA"} = 0.03;  $RSCU_TABLE{"CAG"} = 1.97;
    $RSCU_TABLE{"AAT"} = 0.13;  $RSCU_TABLE{"AAC"} = 1.87;  $RSCU_TABLE{"AAA"} = 0.06;  $RSCU_TABLE{"AAG"} = 1.94;  
    $RSCU_TABLE{"GAT"} = 0.90;  $RSCU_TABLE{"GAC"} = 1.10;  $RSCU_TABLE{"GAA"} = 0.19;  $RSCU_TABLE{"GAG"} = 1.81;
    $RSCU_TABLE{"TGT"} = 0.07;  $RSCU_TABLE{"TGC"} = 1.93;  $RSCU_TABLE{"TGA"} = 0.00;  $RSCU_TABLE{"TGG"} = 1.00;  
    $RSCU_TABLE{"CGT"} = 2.65;  $RSCU_TABLE{"CGC"} = 3.07;  $RSCU_TABLE{"CGA"} = 0.07;  $RSCU_TABLE{"CGG"} = 0.00;  
    $RSCU_TABLE{"AGT"} = 0.04;  $RSCU_TABLE{"AGC"} = 1.13;  $RSCU_TABLE{"AGA"} = 0.00;  $RSCU_TABLE{"AGG"} = 0.21;  
    $RSCU_TABLE{"GGT"} = 1.34;  $RSCU_TABLE{"GGC"} = 1.66;  $RSCU_TABLE{"GGA"} = 0.99;  $RSCU_TABLE{"GGG"} = 0.00;
  }
  elsif ($swit == 6)  #Bacillus subtilis
  {
    $RSCU_TABLE{"TTT"} = 0.70;  $RSCU_TABLE{"TTC"} = 1.30;  $RSCU_TABLE{"TTA"} = 2.71;  $RSCU_TABLE{"TTG"} = 0.00;  
    $RSCU_TABLE{"CTT"} = 2.13;  $RSCU_TABLE{"CTC"} = 0.00;  $RSCU_TABLE{"CTA"} = 1.16;  $RSCU_TABLE{"CTG"} = 0.00;
    $RSCU_TABLE{"ATT"} = 0.91;  $RSCU_TABLE{"ATC"} = 1.96;  $RSCU_TABLE{"ATA"} = 0.13;  $RSCU_TABLE{"ATG"} = 1.00;  
    $RSCU_TABLE{"GTT"} = 1.88;  $RSCU_TABLE{"GTC"} = 0.25;  $RSCU_TABLE{"GTA"} = 1.38;  $RSCU_TABLE{"GTG"} = 0.50;
    $RSCU_TABLE{"TCT"} = 3.45;  $RSCU_TABLE{"TCC"} = 0.00;  $RSCU_TABLE{"TCA"} = 1.50;  $RSCU_TABLE{"TCG"} = 0.00;  
    $RSCU_TABLE{"CCT"} = 2.29;  $RSCU_TABLE{"CCC"} = 0.00;  $RSCU_TABLE{"CCA"} = 1.14;  $RSCU_TABLE{"CCG"} = 0.57;
    $RSCU_TABLE{"ACT"} = 2.21;  $RSCU_TABLE{"ACC"} = 0.00;  $RSCU_TABLE{"ACA"} = 1.38;  $RSCU_TABLE{"ACG"} = 0.41;  
    $RSCU_TABLE{"GCT"} = 2.94;  $RSCU_TABLE{"GCC"} = 0.08;  $RSCU_TABLE{"GCA"} = 0.60;  $RSCU_TABLE{"GCG"} = 0.38;
    $RSCU_TABLE{"TAT"} = 0.50;  $RSCU_TABLE{"TAC"} = 1.50;  $RSCU_TABLE{"TAA"} = 0.00;  $RSCU_TABLE{"TAG"} = 0.00;  
    $RSCU_TABLE{"CAT"} = 2.00;  $RSCU_TABLE{"CAC"} = 0.00;  $RSCU_TABLE{"CAA"} = 1.71;  $RSCU_TABLE{"CAG"} = 0.29;
    $RSCU_TABLE{"AAT"} = 0.47;  $RSCU_TABLE{"AAC"} = 1.53;  $RSCU_TABLE{"AAA"} = 1.83;  $RSCU_TABLE{"AAG"} = 0.17;  
    $RSCU_TABLE{"GAT"} = 0.53;  $RSCU_TABLE{"GAC"} = 1.47;  $RSCU_TABLE{"GAA"} = 1.40;  $RSCU_TABLE{"GAG"} = 0.60;
    $RSCU_TABLE{"TGT"} = 0.00;  $RSCU_TABLE{"TGC"} = 2.00;  $RSCU_TABLE{"TGA"} = 0.00;  $RSCU_TABLE{"TGG"} = 1.00;  
    $RSCU_TABLE{"CGT"} = 3.11;  $RSCU_TABLE{"CGC"} = 1.78;  $RSCU_TABLE{"CGA"} = 0.00;  $RSCU_TABLE{"CGG"} = 0.00;  
    $RSCU_TABLE{"AGT"} = 0.45;  $RSCU_TABLE{"AGC"} = 0.60;  $RSCU_TABLE{"AGA"} = 1.11;  $RSCU_TABLE{"AGG"} = 0.00;  
    $RSCU_TABLE{"GGT"} = 1.38;  $RSCU_TABLE{"GGC"} = 0.97;  $RSCU_TABLE{"GGA"} = 1.66;  $RSCU_TABLE{"GGG"} = 0.00;  
  }
  elsif ($swit == 7)  #Deinococcus radiodurans 10.1016/j.biosystems.2005.12.003
  {
    $RSCU_TABLE{"TTT"} = 0.27;  $RSCU_TABLE{"TTC"} = 1.73;  $RSCU_TABLE{"TTA"} = 0.00;  $RSCU_TABLE{"TTG"} = 0.10;  
    $RSCU_TABLE{"CTT"} = 0.13;  $RSCU_TABLE{"CTC"} = 2.27;  $RSCU_TABLE{"CTA"} = 0.01;  $RSCU_TABLE{"CTG"} = 3.47;
    $RSCU_TABLE{"ATT"} = 0.53;  $RSCU_TABLE{"ATC"} = 2.47;  $RSCU_TABLE{"ATA"} = 0.00;  $RSCU_TABLE{"ATG"} = 1.00;
    $RSCU_TABLE{"GTT"} = 0.11;  $RSCU_TABLE{"GTC"} = 1.39;  $RSCU_TABLE{"GTA"} = 0.04;  $RSCU_TABLE{"GTG"} = 2.46;
    $RSCU_TABLE{"TCT"} = 0.09;  $RSCU_TABLE{"TCC"} = 0.92;  $RSCU_TABLE{"TCA"} = 0.04;  $RSCU_TABLE{"TCG"} = 1.28;
    $RSCU_TABLE{"CCT"} = 0.25;  $RSCU_TABLE{"CCC"} = 2.56;  $RSCU_TABLE{"CCA"} = 0.03;  $RSCU_TABLE{"CCG"} = 1.15;
    $RSCU_TABLE{"ACT"} = 0.14;  $RSCU_TABLE{"ACC"} = 3.26;  $RSCU_TABLE{"ACA"} = 0.03;  $RSCU_TABLE{"ACG"} = 0.58;      
    $RSCU_TABLE{"GCT"} = 0.23;  $RSCU_TABLE{"GCC"} = 2.54;  $RSCU_TABLE{"GCA"} = 0.08;  $RSCU_TABLE{"GCG"} = 1.15;
    $RSCU_TABLE{"TAT"} = 0.11;  $RSCU_TABLE{"TAC"} = 1.89;  $RSCU_TABLE{"TAA"} = 0.00;  $RSCU_TABLE{"TAG"} = 0.00;  
    $RSCU_TABLE{"CAT"} = 0.09;  $RSCU_TABLE{"CAC"} = 1.91;  $RSCU_TABLE{"CAA"} = 0.17;  $RSCU_TABLE{"CAG"} = 1.83;
    $RSCU_TABLE{"AAT"} = 0.08;  $RSCU_TABLE{"AAC"} = 1.92;  $RSCU_TABLE{"AAA"} = 0.22;  $RSCU_TABLE{"AAG"} = 1.78;
    $RSCU_TABLE{"GAT"} = 0.11;  $RSCU_TABLE{"GAC"} = 1.89;  $RSCU_TABLE{"GAA"} = 1.18;  $RSCU_TABLE{"GAG"} = 0.82;
    $RSCU_TABLE{"TGT"} = 0.16;  $RSCU_TABLE{"TGC"} = 1.84;  $RSCU_TABLE{"TGA"} = 0.00;  $RSCU_TABLE{"TGG"} = 1.00;
    $RSCU_TABLE{"CGT"} = 0.84;  $RSCU_TABLE{"CGC"} = 4.67;  $RSCU_TABLE{"CGA"} = 0.04;  $RSCU_TABLE{"CGG"} = 0.38;  
    $RSCU_TABLE{"AGT"} = 0.14;  $RSCU_TABLE{"AGC"} = 3.54;  $RSCU_TABLE{"AGA"} = 0.03;  $RSCU_TABLE{"AGG"} = 0.05;  
    $RSCU_TABLE{"GGT"} = 0.40;  $RSCU_TABLE{"GGC"} = 3.21;  $RSCU_TABLE{"GGA"} = 0.07;  $RSCU_TABLE{"GGG"} = 0.32;  
  }
  elsif ($swit == 8)  #Mycoplasma genitalium 10.1111/j.1472-765X.2004.01619.x, bad numbers for Y and L!!! W (TGG and TGA) from Nakamura:1997p22217
  {
    $RSCU_TABLE{"TTT"} = 1.73;  $RSCU_TABLE{"TTC"} = 0.27;  $RSCU_TABLE{"TTA"} = 2.83;  $RSCU_TABLE{"TTG"} = 0.81;  
    $RSCU_TABLE{"CTT"} = 1.12;  $RSCU_TABLE{"CTC"} = 0.28;  $RSCU_TABLE{"CTA"} = 0.71;  $RSCU_TABLE{"CTG"} = 0.25;
    $RSCU_TABLE{"ATT"} = 1.89;  $RSCU_TABLE{"ATC"} = 0.65;  $RSCU_TABLE{"ATA"} = 0.46;  $RSCU_TABLE{"ATG"} = 1.00;
    $RSCU_TABLE{"GTT"} = 2.47;  $RSCU_TABLE{"GTC"} = 0.22;  $RSCU_TABLE{"GTA"} = 0.85;  $RSCU_TABLE{"GTG"} = 0.46;
    $RSCU_TABLE{"TCT"} = 1.13;  $RSCU_TABLE{"TCC"} = 0.37;  $RSCU_TABLE{"TCA"} = 1.47;  $RSCU_TABLE{"TCG"} = 0.11;
    $RSCU_TABLE{"CCT"} = 1.96;  $RSCU_TABLE{"CCC"} = 0.48;  $RSCU_TABLE{"CCA"} = 1.43;  $RSCU_TABLE{"CCG"} = 0.13;
    $RSCU_TABLE{"ACT"} = 1.89;  $RSCU_TABLE{"ACC"} = 0.76;  $RSCU_TABLE{"ACA"} = 1.23;  $RSCU_TABLE{"ACG"} = 0.12;
    $RSCU_TABLE{"GCT"} = 1.97;  $RSCU_TABLE{"GCC"} = 0.30;  $RSCU_TABLE{"GCA"} = 1.54;  $RSCU_TABLE{"GCG"} = 0.19;
    $RSCU_TABLE{"TAT"} = 1.48;  $RSCU_TABLE{"TAC"} = 0.52;  $RSCU_TABLE{"TAA"} = 0.00;  $RSCU_TABLE{"TAG"} = 0.00;
    $RSCU_TABLE{"CAT"} = 1.30;  $RSCU_TABLE{"CAC"} = 0.70;  $RSCU_TABLE{"CAA"} = 1.62;  $RSCU_TABLE{"CAG"} = 0.38;
    $RSCU_TABLE{"AAT"} = 1.23;  $RSCU_TABLE{"AAC"} = 0.77;  $RSCU_TABLE{"AAA"} = 1.48;  $RSCU_TABLE{"AAG"} = 0.52;
    $RSCU_TABLE{"GAT"} = 1.72;  $RSCU_TABLE{"GAC"} = 0.28;  $RSCU_TABLE{"GAA"} = 1.61;  $RSCU_TABLE{"GAG"} = 0.39;
    $RSCU_TABLE{"TGT"} = 1.59;  $RSCU_TABLE{"TGC"} = 0.41;  $RSCU_TABLE{"TGA"} = 1.30;  $RSCU_TABLE{"TGG"} = 0.70;  
    $RSCU_TABLE{"CGT"} = 1.35;  $RSCU_TABLE{"CGC"} = 0.59;  $RSCU_TABLE{"CGA"} = 0.26;  $RSCU_TABLE{"CGG"} = 0.20;  
    $RSCU_TABLE{"AGT"} = 2.33;  $RSCU_TABLE{"AGC"} = 0.60;  $RSCU_TABLE{"AGA"} = 2.70;  $RSCU_TABLE{"AGG"} = 0.90;  
    $RSCU_TABLE{"GGT"} = 1.99;  $RSCU_TABLE{"GGC"} = 0.44;  $RSCU_TABLE{"GGA"} = 0.99;  $RSCU_TABLE{"GGG"} = 0.58;  
  }
  return \%RSCU_TABLE;
}

sub RSCU_filter
{
  my ($RSCU_TABLE, $min_value) = @_;
  my @bad_codons = grep { $$RSCU_TABLE{$_}  < $min_value } 
           keys %$RSCU_TABLE;
  delete @$RSCU_TABLE{@bad_codons};
  return $RSCU_TABLE;
}

sub define_codon_percentages
{
  my ($CODON_TABLE, $RSCU_VALUES) = @_;
  my %AA_cod_count;
  $AA_cod_count{$$CODON_TABLE{$_}}++  foreach keys %$CODON_TABLE;
  my %CODON_PERC_TABLE = map { $_ => $$RSCU_VALUES{$_} / $AA_cod_count{$$CODON_TABLE{$_}} } keys %$CODON_TABLE; 
  return \%CODON_PERC_TABLE;
}

sub index_codon_percentages
{
  my ($ntseq, $window, $cpthashref) = @_;
  my @xvalues; my @yvalues;
  my %CODON_PERCENTAGE_TABLE = %$cpthashref;
  my $index; my $sum;
  for (my $x = int($window *(3/2))-3; $x < (length($ntseq) - 3*(int($window *(3/2))-3)); $x+=3)
  {
    $sum = 0;
    for(my $y = $x; $y < 3*$window + $x; $y += 3)
    {
      $sum += $CODON_PERCENTAGE_TABLE{substr($ntseq, $y, 3)};
  #    $sum += $RSCU_TABLE{substr($nucseq, $y, 3)};
    }
    $sum = $sum / $window;
    $index = ($x / 3) + 1;
    push @xvalues, $index;
    push @yvalues, $sum;
  }
  return (\@xvalues, \@yvalues);
}

sub codon_count
{
  my ($arrayref, $CODON_TABLE) = @_;
  my %codoncount = map {$_ => 0} keys %$CODON_TABLE;
  foreach my $seq (@$arrayref)
  {
    my $offset = 0;
    while ( $offset <= length($seq)-3 )
    {
      my $codon = substr($seq, $offset, 3);
      if ($codon =~ $strcodon)
      {
        $codoncount{$codon} ++;
      }
      else
      {
        $codoncount{"XXX"} ++;
      }
      $offset += 3;
    }
  }
  return \%codoncount;
}

sub generate_RSCU_values
{
  my ($codon_count, $CODON_TABLE) = @_;
  my $REV_CODON_TABLE = define_reverse_codon_table($CODON_TABLE);
  my $RSCU_hash = {}; 
  foreach (sort grep {$_ ne "XXX"} keys %$codon_count)
  {
    my $x_j = 0;
    my $x = $$codon_count{$_};
    my $family = $$REV_CODON_TABLE{$$CODON_TABLE{$_}};
    my $family_size = scalar(@$family);
    $x_j += $$codon_count{$_} foreach (grep {exists $$codon_count{$_}} @$family);
    $$RSCU_hash{$_} = sprintf("%.3f",  $x / ($x_j / $family_size) ) ;#+ 0;
  }
  return $RSCU_hash;
}

sub define_aa_defaults
{
  my ($CODON_TABLE, $RSCU_VALUES) = @_;
  my $REV_CODON_TABLE = define_reverse_codon_table($CODON_TABLE);
  my %aa_defaults = ();
  foreach my $aa (keys %AACIDS)
  {
    my $myrscu = 0;
    foreach my $codon (@{$$REV_CODON_TABLE{$aa}})  
    {  
      if ($$RSCU_VALUES{$codon} > $myrscu)  
      {
        $aa_defaults{$aa} = $codon;
        $myrscu = $$RSCU_VALUES{$codon};
      }
    }
  }
  return %aa_defaults;
}

1;

__END__

=head1 NAME

GeneDesign::Basic - perl functions for computer assisted synthetic design

=head1 VERSION

Version 3.00

=head1 DESCRIPTION

  GeneDesign is a library for the computer-assisted design of synthetic genes
  
=head1 Functions

=head2 define_aa_defaults()
  Generates a hash.  KEYS: one letter amino acid code (string)  
  VALUES: most highly expressed codon for that amino acid (string)
  in: reverse codon table (hash reference),
      RSCU value table (hash reference)
  out: amino acid default table (hash)

=head2 define_codon_table()
  Generates a hash.  KEYS: codons (string)  VALUES: amino acids (string)
  in: switch for codon table (1)
  out: codon table (hash)

=head2 define_reverse_codon_table()
  Generates a reference to a hash. 
  KEYS: amino acids (string)  VALUES: codon list (array reference)
  in: codon table (hash reference)
  out: reverse codon table (hash)

=head2 define_RSCU_values()
  Generates a hash.  KEYS: codons (string)  VALUES: RSCU value (float)
  in: switch for species (1 s.cer, 2 e.col, 3 h.sap, 4 c.ele, 5 d.mel)
  out: RSCU value table (hash)

=head2 RSCU_filter()
  Deletes anything from the RSCU_TABLE that doesn't meet minimum RSCU.
  in: rscu table (hash reference),
      minimum RSCU value (integer)
  out: rscu table (hash reference)

=head2 define_codon_percentages()
  Generates a hash.  KEYS: codons (string) 
  VALUES: RSCU value over codon family size (float)
  in: codon table (hash reference),
      RSCU value table (hash reference)
  out: codon percentage table (hash)

=head2 index_codon_percentages()
  Generates two arrays for x and y values of a graph of codon percentage values.
  in: dna sequence (string),
      window size (integer),
      codon percentage table (hash reference)
  out: x values (array reference), y values (array reference)

=head2 pattern_remover()
  takes a nucleotide sequence, a nucleotide "pattern" to be removed, and a few 
  codon tables, and returns an edited nucleotide sequence that is missing the 
  pattern (if possible).  Ranks codon replacements by RSCU differences to 
  minimize expression damage.
  in: nucleotide sequence (string), 
      nucleotide pattern (string), 
      codon table (hash reference), 
      RSCU value table (hash reference)
  out: nucleotide sequence (string) OR null

=head2 pattern_adder()
  takes a nucleotide sequence, a nucleotide "pattern" to be interpolated, and 
  the codon table, and returns an edited nucleotide sequence that contains the 
  pattern (if possible).
  in: nucleotide sequence (string), 
      nucleotide pattern (string), 
      codon table (hash reference), 
      RSCU value table (hash reference)
  out: nucleotide sequence (string) OR null

=head2 pattern_aligner()
  takes a nucleotide sequence, a pattern, a peptide sequence, and a codon table
  and inserts Ns before the pattern until they align properly. This is so a
  pattern can be inserted out of frame.
  in: nucleotide sequence (string), 
      nucleotide pattern (string), 
      amino acid sequence (string), 
      codon table (hash reference)
  out: nucleotide pattern (string)

=head2 pattern_finder()

=head2 change_codons()
  takes a nucleotide sequence and a few codon tables and tries to recode the
  nucleotide sequence to one of four algorithms, 0 random, 1 most optimal, 
  2 next most optimal, 3 most different, 4 least different.
  in: nucleotide sequence (string), 
      codon table (hash reference), 
      reverse codon table (hash reference), 
      RSCU value table (hash reference), 
      algorithm number 
      tag: try not to change the first two bases of the first codon in.
  out: nucleotide sequence (string)

=head2 reverse_translate()
  takes an amino acid sequence and a specific codon table and returns that frame
  translated into amino acids.  See gdRevTrans.cgi for use.
  in: nucleotide sequence (string), 
      switch for frame (1, 2, or 3),
      codon table (hash reference)
  out: amino acid sequence (string)

=head2 amb_transcription()
  takes an ambiguous nucleotide sequence and returns a list of all possible non
  ambiguous nucleotide sequences it could represent
  in: nucleotide sequence (string),
      codon table (hash reference),
      reverse codon table (hash reference)
  out: nucleotide sequence list (vector)

=head2 amb_translation()
  takes a nucleotide that may be degenerate and a codon table and returns a list
  of all amino acid sequences that nucleotide sequence could be translated into.
  in: nucleotide sequence (string),
      codon table (hash reference)
  out: amino acid sequence list (vector)

=head2 degcodon_to_aas()
  takes a codon that may be degenerate and a codon table and returns a list of 
  all amino acids that codon could represent. If a hashref is provided with 
  previous answers, it will run MUCH faster.
  in: codon (string), 
      codon table (hash reference)
  out: amino acid list (vector)

=head2 translate()
  takes a nucleotide sequence, a frame, and a codon table and returns that frame
  translated into amino acids.
  in: nucleotide sequence (string),
      switch for frame (±1, ±2, or ±3),
      codon table (hash reference)
  out: amino acid sequence (string)

=head2 codon_count()
  takes a reference to an array of sequences and returns a hash with codons as 
  keys and the number of times the codon occurs as a value.
  in: gene sequence (array reference)
  out: codon count (hash reference)

=head2 generate_RSCU_values()
  takes a hash reference with keys as codons and values as number of times 
  those codons occur (it helps to use codon_count) and returns a hash with each 
  codon and its RSCU value
  in: codon count (hash reference), 
      reverse codon table (hash reference)
  out: RSCU values (hash reference)

=head2 rscu_parser()
  takes the form AAA (K) 0.540

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
