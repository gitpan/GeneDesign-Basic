package GeneDesign::Basic;
require Exporter;

use POSIX qw(log10);
use List::Util qw(shuffle first);
use Class::Struct;
use Text::Wrap qw($columns &wrap);

use strict;
use warnings;

use vars qw($VERSION @ISA @EXPORT_OK %EXPORT_TAGS);
$VERSION = 3.00;

@ISA = qw(Exporter);
@EXPORT_OK = qw(
  regres 
  make_oligos
  compare_sequences
  count
  ntherm
  compareseqs
  regres
  complement
  melt
  cleanup
  oligocruncher
  orf_finder
  define_oligos
  fasta_parser
  cons_seq
  print_alignment
  fasta_writer
  %NTIDES
  %AA_NAMES
  $ambnt
  %ORGANISMS
  $docpath
  $linkpath
  $enzfile
  %AACIDS
);

%EXPORT_TAGS =  (all => [qw(regres make_oligos compare_sequences
  count ntherm compareseqs translate regres complement melt cleanup 
  oligocruncher orf_finder define_oligos fasta_parser cons_seq print_alignment
  fasta_writer %NTIDES %AA_NAMES $ambnt %ORGANISMS %AACIDS $docpath $linkpath 
  $enzfile)]);

struct Oligo => { map { $_ => '$' } qw(ChunkNumber OligoNumber OligoLength OligoStart OligoStop OligoSense FpOlapLeng TpOlapLeng OligoSeq) };

struct Chunk => { map { $_ => '$' } qw(Parameters ChunkNumber NumberofOligos ChunkLength AvgOligoLength ChunkStart ChunkStop Collisions AvgGapLeng AvgOlapLeng 
          AvgOlapMelt ChunkSeq Oligos Olaps Users FivePrimeEnz ThreePrimeEnz ShrtOligo LongOligo ThreePrimeOlap Mask Name) };

struct USERsite => { map { $_ => '$' } qw(Start nNumber Sequence) };

our $docpath = "../../Documents/gd";
our $linkpath = "http://localhost/gd";
our $enzfile = "bs_enzymes.txt";

our %NTIDES = (A => "A", B => "[BCGKSTY]", C => "C", D => "[ADGKRTW]", G => "G", H => "[ACHMTWY]", K => "[GKT]", M => "[ACM]", 
           N => "[ABCDGHKMNRSTVWY]", R => "[AGR]", S => "[CGS]", T => "T", V => "[ACGMRSV]", W => "[ATW]", Y => "[CTY]", );
our @NTS = qw(A T C G);
our @NTSG = qw(B D H K M N R S V W Y);

our %AACIDS = map { $_, $_ } qw(A C D E F G H I K L M N P Q R S T V W Y);
$AACIDS{"*"} = "[\*]";

our %AA_NAMES = (A => "Ala", B => "Unk", C => "Cys", D => "Asp", E => "Glu", F => "Phe", G => "Gly", H => "His", I => "Ile",  
         J => "Unk", K => "Lys", L => "Leu", M => "Met", N => "Asn", O => "Unk", P => "Pro", Q => "Gln", R => "Arg", 
         S => "Ser", T => "Thr", U => "Unk", V => "Val", W => "Trp", X => "Unk", Y => "Tyr", Z => "Unk","*" => "Stp");

# entropy, enthalpy, and free energy of paired bases ÆHû   ÆSû  ÆGû
my %TEEFE = ("TC" => ([ 8.8, 23.5, 1.5]), "GA" => ([ 8.8, 23.5, 1.5]), "CT" => ([ 6.6, 16.4, 1.5]), "AG" => ([ 6.6, 16.4, 1.5]),
       "GG" => ([10.9, 28.4, 2.1]), "CC" => ([10.9, 28.4, 2.1]), "AA" => ([ 8.0, 21.9, 1.2]), "TT" => ([ 8.0, 21.9, 1.2]),
       "AT" => ([ 5.6, 15.2, 0.9]), "TA" => ([ 6.6, 18.4, 0.9]), "CG" => ([11.8, 29.0, 2.8]), "GC" => ([10.5, 26.4, 2.3]),
       "CA" => ([ 8.2, 21.0, 1.7]), "TG" => ([ 8.2, 21.0, 1.7]), "GT" => ([ 9.4, 25.5, 1.5]), "AC" => ([ 9.4, 25.5, 1.5]));

our %ORGANISMS = (0 => "(no organism defined)", 1 => "Saccharomyces cerevisiae", 2 => "E. coli", 
          3 => "Homo sapiens", 4 => "C. elegans", 5 => "Drosophila melanogaster", 6 => "Bacillus subtilis",
          7 => "Deinococcus radiodurans", 8 => "Mycoplasma genitalium");
  
my $nonawords  = qr/[^\w*]*/;
my $nonwords  = qr/\W*/;
my $fastaline  = qr/\A>[\S\ ]*[\n\r]{1}/;
my $nonaa    = qr/[BJOUX]{1}/;
our $ambnt    = qr/[RYWSKMBDHVN]+/;

################################################################################
###########################                          ###########################
################################################################################

sub count
{
  my ($strand) = @_;
  return if (!$strand);
  my $BC = {};
  $$BC{'length'} = length($strand);
  foreach (@NTS, @NTSG)  {  $$BC{$_}   = ($strand =~ s/$_//ig || 0);  };
  foreach (@NTS)      {  $$BC{'d'} += $$BC{$_};      };
  foreach (@NTSG)      {  $$BC{'n'} += $$BC{$_};      };
  $$BC{'?'} = ($$BC{'d'} + $$BC{'n'}) - $$BC{'length'};
  $$BC{'U'} = ($strand =~ s/U//ig || 0);
  my $split = .5*$$BC{'R'}    + .5*$$BC{'Y'}    + .5*$$BC{'K'}    + .5*$$BC{'M'}    + .5*$$BC{'N'};
  my $trip  = (2/3)*$$BC{'B'} + (2/3)*$$BC{'V'} + (1/3)*$$BC{'D'} + (1/3)*$$BC{'H'};
  if ($trip || $split)
  {
    $$BC{'GCp'} = int(((($$BC{'S'}+$$BC{'G'}+$$BC{'C'}+$split + $trip)/$$BC{'length'})*100)+.5);
    $$BC{'ATp'} = int(((($$BC{'W'}+$$BC{'A'}+$$BC{'T'}+$split + (.9-$trip))/$$BC{'length'})*100)+.5);
  }
  elsif (! $trip && ! $split)
  {
    $$BC{'GCp'} = int(((($$BC{'S'}+$$BC{'G'}+$$BC{'C'})/$$BC{'length'})*100)+.5);
    $$BC{'ATp'} = int(((($$BC{'W'}+$$BC{'A'}+$$BC{'T'})/$$BC{'length'})*100)+.5);
  }
  return $BC;
}

sub melt
{
  my ($strand, $swit, $salt, $conc) = @_;
  return if (!$strand);
  $swit = $swit || 3;
  $salt = $salt || .05;
  $conc = $conc || .0000001;
  my $mgc = 1.987;
  my $BC = count($strand);
  if ($swit == 1) #simple
  {
    return ((4 * ($$BC{'C'} + $$BC{'G'})) + (2 * ($$BC{'A'} + $$BC{'T'})));
  }
  if ($swit == 2 || $swit == 3) #baldwin, primer3
  {
    my $base = 81.5 + 16.6*log10($salt) + 41*(($$BC{'C'}+$$BC{'G'})/length($strand));
    return $base - (675/length($strand)) if ($swit == 2);
    return $base - (600/length($strand)) if ($swit == 3);
  }
  if ($swit == 4) ##nntherm
  {
    my ($dH, $dS, $dG) = ntherm($strand);
    return  ((($dH-3.4) / (($dS+($mgc*abs(log($conc / 2))))/1000))-273.15) + (16.6*log10($salt));
  }
  return undef;
}

sub ntherm
{
  my ($strand) = @_;
  my ($dH, $dS, $dG) = (0, 0, 0);
  foreach my $w (keys %TEEFE)
  {  
    while ($strand =~ /(?=$w)/ig)
    {
      $dH += $TEEFE{$w}->[0];
      $dS += $TEEFE{$w}->[1];
      $dG += $TEEFE{$w}->[2];
    }
  }
  return ($dH, $dS, $dG);
}

sub complement
{
  my ($strand, $swit) = @_;
  return undef if (!$strand);
  $strand = scalar reverse($strand) if ($swit);
  $strand =~ tr/ACGTRYKMBDHV/TGCAYRMKVHDB/;
  return $strand;
}

sub regres
{
  my ($sequence, $swit) = @_;
  return if (!$sequence);
  $swit = 1 if (!$swit);
  my $comp = "";
  foreach my $char (split('', $sequence))
  {
    if ($swit == 1)
    {
      $comp .= exists $NTIDES{$char}  ?  $NTIDES{$char}  :  "[X]";
    }
    elsif ($swit == 2)
    {
      $comp .= exists $AACIDS{$char}  ?  $AACIDS{$char}  :  "[X]";
    }  
  }
  return $comp;
}

sub compareseqs
{
  my ($cur, $tar) = @_;
  return 1 if ($tar =~ regres($cur, 1) || $cur =~ regres($tar, 1));
  return 0;
}

sub cleanup
{
  my ($sequence, $swit) = @_;
  $swit = 0 if (!$swit);
  $sequence =~ s/$fastaline//;            #remove FASTA info line
  $sequence = uc $sequence;              #capitalize everything
  $sequence =~ s/$nonawords//g  if ($swit == 2);  #remove every nonword except *
  $sequence =~ s/$nonwords//g    if ($swit != 2);  #remove every nonword
  if ($swit < 2)  {  $sequence =~ s/$_//g foreach ( qw(E F I J L O P Q U X Z) );  }  #strict nucleotide editing
  if ($swit == 0)  {  $sequence =~ s/$_//g foreach (  @NTSG           );  }  #degenerate nucleotide editing
  if ($swit == 2)  {  $sequence =~ s/$_//g foreach ( qw(B J O U X Z)       );  }  #amino acid editing
  return $sequence;
}

sub compare_sequences
{
  my ($topseq, $botseq) = @_;
  return if (!$botseq || length($botseq) != length($topseq));
  my ($tsit, $tver, $len) = (0, 0, length($topseq));
  my $alresults;
  while (length($topseq) > 0)
  {
    my ($topbit, $botbit) = (chop($topseq), chop ($botseq));
    if ($topbit ne $botbit)
    {
      $topbit = $topbit =~ $NTIDES{R}  ?  1  :  0;
      $botbit = $botbit =~ $NTIDES{R}  ?  1  :  0;
      $tsit++ if ($topbit == $botbit);
      $tver++ if ($topbit != $botbit);
    }
  }
  $$alresults{'D'} = $tsit + $tver;                    #changes
  $$alresults{'I'} = $len - $$alresults{'D'};                #identities
  $$alresults{'T'} = $tsit;                        #transitions
  $$alresults{'V'} = $tver;                        #transversions
  $$alresults{'P'} = int((100 - (($$alresults{'D'} / $len) * 100)) + .5);  #percent identity
  return $alresults;
}

sub cons_seq
{
  my ($arrref) = @_;
  my $cons = '';
  for my $x (0..length($$arrref[0]))
  {
    my $flag = 0;
    my $init = substr($$arrref[0], $x, 1);
    for my $y (1..scalar(@$arrref)-1)
    {
      $flag++ if ($init ne substr($$arrref[$y], $x, 1));
    }
    $cons .=  $flag == 0  ?  substr($$arrref[0], $x, 1)  :  '_';
  }
  return $cons;
}

sub fasta_parser
{
  my ($instr, $swit) = @_;
  $swit = 0 if (!$swit);
  my $seqhsh = {};
  my @pre = split(">", $instr);
  shift @pre;
  foreach my $preseq (@pre)
  {
    my @pair = split(/[\n\r]/, $preseq);
    my $id = shift @pair;
    if ($swit)
    {
      my @preid = split(/\s/, $id);
      $id = $preid[0];
    }
    $$seqhsh{">" . $id} = uc join("", @pair);
  }
  return $seqhsh;
}

sub fasta_writer
{
  my ($seqhsh) = @_;
  my $outstr = '';
  $columns = 81;
  foreach my $id (sort {$a cmp $b} keys %$seqhsh)
  {
    $outstr .= $id . "\n";
    $outstr .= wrap("","", $$seqhsh{$id}). "\n";
  }
  return $outstr;
}

sub print_alignment
{
  my ($nuchshref, $width, $swit, $aaseq) = @_;
  my ($space, $break, $times) = $swit == 0  ?  (" ", "\n", "x")  :  ("&nbsp;", "<br>", "&times;");
  my $start = 0;
  my $output = $break;
  while ($start < length($$nuchshref{first {1} keys %$nuchshref}))
  {
    my $count = 0;
    my $iter = 1;
    if ($aaseq)
    {
      my $rt;
      $output.= $break . "0 tran" . $space;
      for ($rt = ($start / 3); $rt < ($start / 3) + ($width / 3); $rt++)
      {
        $output .= substr($aaseq, $rt, 1) . $space x 2 if ($rt < length($aaseq));
        $output .= $space x 3 if ($rt >= length($aaseq));
      }
      $output .= $space . $rt;
    }
    foreach (sort keys %$nuchshref)
    {
      $output .= $break . $_ . $space . substr($$nuchshref{$_}, $start, $width);
      $output .= $space . ($start+$width)*$iter if ($count == 0 && ($start+$width)*$iter < length($$nuchshref{$_}));
      $output .= ($space x (($start+$width)*$iter - length($$nuchshref{$_})+1)) . ($start+$width)*$iter if ($count == 0 && ($start+$width)*$iter >= length($$nuchshref{$_}));
      $count++;
    }
    $iter++;
    $output .= $break . $space x 5 . (($space x 9 . $times) x ($width / 10));
    $output .= $break x 2;
    $start += $width;
  }
  return $output;
}

sub oligocruncher
{
  my ($tov, $hashref) = @_;
  my ($tar_chn_len, $tar_cur_dif, $cur_oli_num, $cur_oli_lap, $cur_oli_len, $cur_chn_mel, $cur_oli_gap, $avg_chn_mel, $avg_oli_len, $start, $starte, $starto, $avg) = 0;
  my (@Overlaps, @tree, @begs, @ends, @Oligos);
  my %pa = %$hashref;  
  my %Collisions;
  $tar_chn_len = $pa{per_chn_len};  
  $cur_oli_num = $pa{tar_oli_num};  
  $cur_oli_len = $pa{tar_oli_len};
  $cur_oli_gap = $pa{tar_oli_gap};  
  $cur_oli_lap = $pa{tar_oli_lap};  
  $tar_cur_dif = $tov->ChunkLength - $tar_chn_len;
#print "\n<br><br>\nrev 0 ", $tov->ChunkNumber, ", $tar_chn_len bp, dif $tar_cur_dif, num $cur_oli_num, len $cur_oli_len, lap $cur_oli_lap, mel $pa{tar_chn_mel}<br>";
  if (abs($tar_cur_dif) >= ($pa{tar_oli_len}+$pa{tar_oli_gap}))  ##-if difference btw perfect and current is bigger than another pair of oligos, increment oli_num
  {
    $cur_oli_num = $pa{tar_oli_num} + (2 * int(($tar_cur_dif / ($pa{tar_oli_len} + $pa{tar_oli_gap})) + .5 ));
    $tar_chn_len = ($cur_oli_num * ($pa{tar_oli_gap} + $pa{tar_oli_lap})) + $pa{tar_oli_lap};
    $tar_cur_dif = $tov->ChunkLength - $tar_chn_len;
  }
  $tov->ShrtOligo(2*$pa{max_oli_len});  $tov->LongOligo(0);  
#print "rev 1 ", $tov->ChunkNumber, ", per $tar_chn_len, dif $tar_cur_dif, num $cur_oli_num, len $cur_oli_len, lap $cur_oli_lap, mel $pa{tar_chn_mel}, tol $pa{chn_mel_tol}, <br>";  
  if ($pa{gapswit} == 1)
  {
    ##-if difference can be spread equally across oligos, increase length
    if (abs($tar_cur_dif) >= $cur_oli_num)          
    {
      $cur_oli_len = $pa{tar_oli_len} + int($tar_cur_dif / $cur_oli_num);
      $tar_cur_dif = $tar_cur_dif - ($cur_oli_num * (int($tar_cur_dif / $cur_oli_num)));  
    }
    ##-if the length is violating max_len, increase num by 2, decrease len by 10, recalc
    if ( ($cur_oli_len >= $pa{max_oli_len}) || ($cur_oli_len == $pa{max_oli_len} && $tar_cur_dif > 0) )    
    {
      $cur_oli_len = $pa{tar_oli_len}-10; 
      $cur_oli_num += 2; 
      $tar_chn_len = $cur_oli_num*($cur_oli_len - $pa{tar_oli_lap}) + $pa{tar_oli_lap};
      $tar_cur_dif = $tov->ChunkLength - $tar_chn_len;
      if (abs($tar_cur_dif) >= $cur_oli_num)
      {
        $cur_oli_len = $cur_oli_len + int($tar_cur_dif / $cur_oli_num);
        $tar_cur_dif = $tar_cur_dif - ($cur_oli_num * (int($tar_cur_dif / $cur_oli_num)));          
      }
    }      
    if ($cur_oli_len >= $pa{max_oli_len} && $tar_cur_dif > 0)
    {
      print "oh no, after rev1, current target is >= the max! $pa{max_oli_len} --- please tell Sarah the following string.<br> ";
      print "&nbsp;&nbsp;cur_oli_len $cur_oli_len, cur_oli_num $cur_oli_num, tar_cur_dif $tar_cur_dif, tar_chn_len $tar_chn_len, chunk length ", $tov->ChunkLength, "<br><br>";
    }
    my $start = 0;      
    for (my $w = 1; $w <= $cur_oli_num; $w++)          ##-difference now be between 0 and abs(oli_num-1) - so in/decrement individual overlap lengths
    {
      my $strlen = $cur_oli_len;
      $strlen++ if ( $w <= abs($tar_cur_dif) && $tar_cur_dif > 0);
      $strlen-- if ( $w <= abs($tar_cur_dif) && $tar_cur_dif < 0);
      push @Overlaps, substr($tov->ChunkSeq, $start + $strlen - $cur_oli_lap, $cur_oli_lap) if ($w != $cur_oli_num);
      $start =  $start + $strlen - $cur_oli_lap;
    }
    @tree = map (melt($_, $pa{melform} , .05, .0000001), @Overlaps);
    foreach (@tree)  {  $avg += $_;  }  $avg = int(($avg / scalar(@tree))+.5);
    $cur_chn_mel = ($avg < ($pa{tar_chn_mel}-10))  ?  $pa{tar_chn_mel} - .2*($pa{tar_chn_mel}-$avg) :  $pa{tar_chn_mel};  #Adjust target melting temp for reality.
#print "rev 3 ", $tov->ChunkNumber, ", per $tar_chn_len, dif $tar_cur_dif, num $cur_oli_num, len $cur_oli_len, lap $cur_oli_lap, mel $cur_chn_mel, tol $pa{chn_mel_tol}, <br>";  

    $start = 0;
    @Overlaps = ();
    for (my $w = 1; $w <= $cur_oli_num; $w++)          ##-then make oligos, changing overlaps for melting temperature if appropriate
    {
      my $laplen = $cur_oli_lap;
      my $strlen = $cur_oli_len;
      $strlen++ if ( $w <= abs($tar_cur_dif) && $tar_cur_dif > 0);
      $strlen-- if ( $w <= abs($tar_cur_dif) && $tar_cur_dif < 0);
      $laplen-- if ($strlen < $pa{tar_oli_len} && $cur_oli_len < $pa{tar_oli_len});#
      if ($w != $cur_oli_num)
      {
        unless ($pa{hardlap})
        {
  #    print ".rev3. ", $tov->ChunkNumber, ", $w strlen $strlen, laplen $laplen, num $cur_oli_num, len $cur_oli_len, mel $cur_chn_mel, tol $pa{chn_mel_tol}, max $pa{max_oli_len}, <br>";
          while (melt(substr($tov->ChunkSeq, $start + $strlen - $laplen, $laplen), $pa{melform}, .05, .0000001) >= ($cur_chn_mel + $pa{chn_mel_tol}) && $strlen > $cur_oli_len)
          {
  #      print "...deccing<br>";
            $laplen--;  $strlen--;
          }
          while (melt(substr($tov->ChunkSeq, $start + $strlen - $laplen, $laplen), $pa{melform}, .05, .0000001) <= ($cur_chn_mel - $pa{chn_mel_tol}) && $strlen < $pa{max_oli_len})
          {
  #      print "...inccing<br>";
            $laplen++;  $strlen++;
          }
        }
        push @Overlaps, substr($tov->ChunkSeq, $start + $strlen - $laplen, $laplen);
#      print "&nbsp;&nbsp;$Overlaps[-1]<Br>";
        $avg_chn_mel += melt(substr($tov->ChunkSeq, $start + $strlen - $laplen, $laplen), $pa{melform}, .05, .0000001);
      }
      push @Oligos, substr($tov->ChunkSeq, $start, $strlen);
      push @begs, $start;
      push @ends, $start + length($Oligos[-1]);
      $tov->ShrtOligo(length($Oligos[-1])) if length($Oligos[-1]) <= $tov->ShrtOligo;
      $tov->LongOligo(length($Oligos[-1])) if length($Oligos[-1]) >= $tov->LongOligo;
      $avg_oli_len += length($Oligos[-1]);
      $start =  $start + $strlen - $laplen;
    }
    for (my $w = 0; $w <= $cur_oli_num - 3; $w++)  {  $Collisions{$w} = $ends[$w] - $begs[$w+2]  if ($ends[$w] > $begs[$w+2]);  }
  }
  elsif ($pa{gapswit} == 0)
  {
    $cur_oli_len = int((2* $tov->ChunkLength) / ((2*$tar_chn_len) / $pa{tar_oli_len}));
    $cur_oli_lap = int($cur_oli_len * .5);
    $tar_cur_dif = $tov->ChunkLength - ($cur_oli_num * $cur_oli_len * .5 + $cur_oli_lap);
  
    $starto = $cur_oli_lap + 1;
    for (my $w = 1; $w <= $cur_oli_num; $w++)          ##-difference should now be between 0 and abs(oli_num-1) - so in/decrement individual overlap lengths
    {      
      my $strlen = $cur_oli_len;
      if ( $w <= abs(2 * $tar_cur_dif) && $tar_cur_dif > 0)  {$strlen++;}
      if ( $w <= abs(2 * $tar_cur_dif) && $tar_cur_dif < 0)  {$strlen--;}
      if ($w % 2 == 1)
      {
        push @Oligos, substr($tov->ChunkSeq, $starte, $strlen);
        $starte = $starte + $strlen;
        push @Overlaps, substr($tov->ChunkSeq, $starto, $starte - $starto);
        $avg_chn_mel += melt(substr($tov->ChunkSeq, $starto, $starte - $starto), $pa{melform}, .05, .0000001);
      }
      if ($w % 2 == 0)
      {
        push @Oligos, substr($tov->ChunkSeq, $starto, $strlen);
        $starto = $starto + $strlen;
        push @Overlaps, substr($tov->ChunkSeq, $starte, $starto - $starte);
        $avg_chn_mel += melt(substr($tov->ChunkSeq, $starte, $starto - $starte), $pa{melform}, .05, .0000001);
      }
      $tov->ShrtOligo(length($Oligos[-1])) if length($Oligos[-1]) <= $tov->ShrtOligo;
      $tov->LongOligo(length($Oligos[-1])) if length($Oligos[-1]) >= $tov->LongOligo;
      $avg_oli_len += length($Oligos[-1]);
    }
  }
  $tov->Collisions(\%Collisions);  
  $tov->AvgOlapMelt(int(($avg_chn_mel / scalar(@Overlaps))+.5));
  $tov->AvgOligoLength(int(($avg_oli_len / scalar(@Oligos))+.5));
  $tov->Oligos(\@Oligos);
  $tov->Olaps(\@Overlaps);
  return;
}

sub make_oligos
{
  my ($bbseq, $pa) = @_;
  my ($tar_chn_len, $tar_cur_dif, $cur_oli_num, $cur_oli_lap, $cur_oli_len, 
    $cur_chn_mel, $cur_oli_gap, $avg_chn_mel, $avg_oli_len, 
    $start, $starte, $starto, $avg) = 0;
  my (@Overlaps, @tree, @begs, @ends, @Oligos);
  my %Collisions;
  my $bblen = length($bbseq);
  my $bbinfo = {};
  my $loxp = "ATAACTTCGTATAATGTACATTATACGAAGTTAT";
  
  $tar_chn_len = $$pa{per_chn_len};  
  $cur_oli_num = $$pa{tar_oli_num};  
  $cur_oli_len = $$pa{tar_oli_len};
  $cur_oli_gap = $$pa{tar_oli_gap};  
  $cur_oli_lap = $$pa{tar_oli_lap};  
  $tar_cur_dif = $bblen - $tar_chn_len;
#print "\n<br><br>\nrev 0 ", $tov->ChunkNumber, ", $tar_chn_len bp, dif $tar_cur_dif, num $cur_oli_num, len $cur_oli_len, lap $cur_oli_lap, mel $$pa{tar_chn_mel}<br>";
  if (abs($tar_cur_dif) >= ($$pa{tar_oli_len}+$$pa{tar_oli_gap}))  ##-if difference btw perfect and current is bigger than another pair of oligos, increment oli_num
  {
    $cur_oli_num = $$pa{tar_oli_num} + (2 * int(($tar_cur_dif / ($$pa{tar_oli_len} + $$pa{tar_oli_gap})) + .5 ));
    $tar_chn_len = ($cur_oli_num * ($$pa{tar_oli_gap} + $$pa{tar_oli_lap})) + $$pa{tar_oli_lap};
    $tar_cur_dif = $bblen - $tar_chn_len;
  }
  $bbinfo->{SHORTOLIGO} = 2*$$pa{max_oli_len};
  $bbinfo->{LONGOLIGO} = 0;
#print "rev 1 ", $tov->ChunkNumber, ", per $tar_chn_len, dif $tar_cur_dif, num $cur_oli_num, len $cur_oli_len, lap $cur_oli_lap, mel $$pa{tar_chn_mel}, tol $$pa{chn_mel_tol}, <br>";  
  if ($$pa{gapswit} == 1)
  {
    ##-if difference can be spread equally across oligos, increase length
    if (abs($tar_cur_dif) >= $cur_oli_num)          
    {
      $cur_oli_len = $$pa{tar_oli_len} + int($tar_cur_dif / $cur_oli_num);
      $tar_cur_dif = $tar_cur_dif - ($cur_oli_num * (int($tar_cur_dif / $cur_oli_num)));  
    }
    ##-if the length is violating max_len, increase num by 2, decrease len by 10, recalc
    if ( ($cur_oli_len >= $$pa{max_oli_len}) || ($cur_oli_len == $$pa{max_oli_len} && $tar_cur_dif > 0) )    
    {
      $cur_oli_len = $$pa{tar_oli_len}-10; 
      $cur_oli_num += 2; 
      $tar_chn_len = $cur_oli_num*($cur_oli_len - $$pa{tar_oli_lap}) + $$pa{tar_oli_lap};
      $tar_cur_dif = $bblen - $tar_chn_len;
      if (abs($tar_cur_dif) >= $cur_oli_num)
      {
        $cur_oli_len = $cur_oli_len + int($tar_cur_dif / $cur_oli_num);
        $tar_cur_dif = $tar_cur_dif - ($cur_oli_num * (int($tar_cur_dif / $cur_oli_num)));          
      }
    }      
    if ($cur_oli_len >= $$pa{max_oli_len} && $tar_cur_dif > 0)
    {
      print "oh no, after rev1, current target is >= the max! $$pa{max_oli_len} --- please tell Sarah the following string.<br> ";
      print "&nbsp;&nbsp;cur_oli_len $cur_oli_len, cur_oli_num $cur_oli_num, tar_cur_dif $tar_cur_dif, tar_chn_len $tar_chn_len, chunk length ", $bblen, "<br><br>";
    }
    my $start = 0;      
    for (my $w = 1; $w <= $cur_oli_num; $w++)          ##-difference now be between 0 and abs(oli_num-1) - so in/decrement individual overlap lengths
    {
      my $strlen = $cur_oli_len;
      $strlen++ if ( $w <= abs($tar_cur_dif) && $tar_cur_dif > 0);
      $strlen-- if ( $w <= abs($tar_cur_dif) && $tar_cur_dif < 0);
      push @Overlaps, substr($bbseq, $start + $strlen - $cur_oli_lap, $cur_oli_lap) if ($w != $cur_oli_num);
      $start =  $start + $strlen - $cur_oli_lap;
    }
    @tree = map (melt($_, $$pa{melform} , .05, .0000001), @Overlaps);
    foreach (@tree)  {  $avg += $_;  }  $avg = int(($avg / scalar(@tree))+.5);
    $cur_chn_mel = ($avg < ($$pa{tar_chn_mel}-10))  ?  $$pa{tar_chn_mel} - .2*($$pa{tar_chn_mel}-$avg) :  $$pa{tar_chn_mel};  #Adjust target melting temp for reality.
#print "rev 3 ", $tov->ChunkNumber, ", per $tar_chn_len, dif $tar_cur_dif, num $cur_oli_num, len $cur_oli_len, lap $cur_oli_lap, mel $cur_chn_mel, tol $$pa{chn_mel_tol}, <br>";  

    $start = 0;
    @Overlaps = ();
    for (my $w = 1; $w <= $cur_oli_num; $w++)          ##-then make oligos, changing overlaps for melting temperature if appropriate
    {
      my $laplen = $cur_oli_lap;
      my $strlen = $cur_oli_len;
      $strlen++ if ( $w <= abs($tar_cur_dif) && $tar_cur_dif > 0);
      $strlen-- if ( $w <= abs($tar_cur_dif) && $tar_cur_dif < 0);
      $laplen-- if ($strlen < $$pa{tar_oli_len} && $cur_oli_len < $$pa{tar_oli_len});#
      if ($w != $cur_oli_num)
      {
      
        unless ($$pa{hardlap})
        {
  #    print ".rev3.  $w strlen $strlen, laplen $laplen, num $cur_oli_num, len $cur_oli_len, mel $cur_chn_mel, tol $$pa{chn_mel_tol}, max $$pa{max_oli_len}\n";
          while (melt(substr($bbseq, $start + $strlen - $laplen, $laplen), $$pa{melform}, .05, .0000001) 
            >= ($cur_chn_mel + $$pa{chn_mel_tol}) && $strlen > $cur_oli_len)
          {
  #      print "...deccing<br>";
            $laplen--;  $strlen--;
          }
          while (melt(substr($bbseq, $start + $strlen - $laplen, $laplen), $$pa{melform}, .05, .0000001) 
            <= ($cur_chn_mel - $$pa{chn_mel_tol}) && $strlen < $$pa{max_oli_len})
          {
  #      print "...inccing<br>";
            $laplen++;  $strlen++;
          }   
           if (substr($bbseq, $start, $strlen) =~ /\A$loxp/)
            {
              $start += 5; $strlen -= 5;
            }
             elsif (substr($bbseq, $start, $strlen) =~ /$loxp\Z/)
              {
                $start -= 5; $strlen -= 5;
              }
        }
        push @Overlaps, substr($bbseq, $start + $strlen - $laplen, $laplen);
#      print "\t\t$Overlaps[-1]\n";
        $avg_chn_mel += melt(substr($bbseq, $start + $strlen - $laplen, $laplen), $$pa{melform}, .05, .0000001);
      }
      push @Oligos, substr($bbseq, $start, $strlen);
      push @begs, $start;
      push @ends, $start + length($Oligos[-1]);
      $bbinfo->{SHORTOLIGO} = length($Oligos[-1]) if length($Oligos[-1]) < $bbinfo->{SHORTOLIGO};
      $bbinfo->{LONGOLIGO} = length($Oligos[-1]) if length($Oligos[-1]) > $bbinfo->{LONGOLIGO};
      $avg_oli_len += length($Oligos[-1]);
      $start =  $start + $strlen - $laplen;
    }
    for (my $w = 0; $w <= $cur_oli_num - 3; $w++)  {  $Collisions{$w} = $ends[$w] - $begs[$w+2]  if ($ends[$w] > $begs[$w+2]);  }
  }
  elsif ($$pa{gapswit} == 0)
  {
    $cur_oli_len = int((2* $bblen) / ((2*$tar_chn_len) / $$pa{tar_oli_len}));
    $cur_oli_lap = int($cur_oli_len * .5);
    $tar_cur_dif = $bblen - ($cur_oli_num * $cur_oli_len * .5 + $cur_oli_lap);
  
    $starto = $cur_oli_lap + 1;
    for (my $w = 1; $w <= $cur_oli_num; $w++)          ##-difference should now be between 0 and abs(oli_num-1) - so in/decrement individual overlap lengths
    {      
      my $strlen = $cur_oli_len;
      if ( $w <= abs(2 * $tar_cur_dif) && $tar_cur_dif > 0)  {$strlen++;}
      if ( $w <= abs(2 * $tar_cur_dif) && $tar_cur_dif < 0)  {$strlen--;}
      if ($w % 2 == 1)
      {
        push @Oligos, substr($bbseq, $starte, $strlen);
        $starte = $starte + $strlen;
        push @Overlaps, substr($bbseq, $starto, $starte - $starto);
        $avg_chn_mel += melt(substr($bbseq, $starto, $starte - $starto), $$pa{melform}, .05, .0000001);
      }
      if ($w % 2 == 0)
      {
        push @Oligos, substr($bbseq, $starto, $strlen);
        $starto = $starto + $strlen;
        push @Overlaps, substr($bbseq, $starte, $starto - $starte);
        $avg_chn_mel += melt(substr($bbseq, $starte, $starto - $starte), $$pa{melform}, .05, .0000001);
      }
      $bbinfo->{SHORTOLIGO} = length($Oligos[-1]) if length($Oligos[-1]) < $bbinfo->{SHORTOLIGO};
      $bbinfo->{LONGOLIGO} = length($Oligos[-1]) if length($Oligos[-1]) > $bbinfo->{LONGOLIGO};
      $avg_oli_len += length($Oligos[-1]);
    }
  }
  $bbinfo->{COLLISIONS} = \%Collisions;
  $bbinfo->{AVGOLAPMELT} = int(($avg_chn_mel / scalar(@Overlaps))+.5);
  $bbinfo->{AVGOLIGOLEN} = int(($avg_oli_len / scalar(@Oligos))+.5);
  $bbinfo->{OLIGOS} = \@Oligos;
  $bbinfo->{OLAPS} = \@Overlaps;
  return $bbinfo;
}

sub orf_finder
{
  my ($strand, $hashref) = @_;
  my $answer = [];
  for my $frame (qw(1 2 3 -1 -2 -3))
  {
    my $strandaa = translate($strand, $frame, $hashref);
    my $leng = length($strandaa);
    my $curpos = 0; 
    my $orflength = 0; 
    my $onnaorf = 0; 
    my $orfstart = 0; 
    while ($curpos <= $leng)
    {
      my $aa = substr($strandaa, $curpos, 1);
      if ($aa eq 'M' && $onnaorf eq '0')
      {
        $onnaorf = 1;
        $orfstart = $curpos;
      }
      if ($aa eq '*' || ($curpos == $leng && $onnaorf == 1))
      {
        $onnaorf= 0;
        push @$answer, [$frame, $orfstart, $orflength] if ($orflength >= .1*($leng));
        $orflength = 0;
      }
      $curpos++;
      $orflength++ if ($onnaorf == 1);
    }
  }
  return $answer;
}

sub define_oligos
{
  my ($ollist_ref, $revcompswit) = @_;
  my %OL_DATA;
  ( $OL_DATA{REGEX},  $OL_DATA{CLEAN} ) = ( {}, {} );
  foreach my $oligo (@$ollist_ref)
  {
    $OL_DATA{CLEAN}->{$oligo}  = $oligo;    #recognition site
  #regular expression array  
    my @arr = ( $revcompswit == 1 && complement($oligo, 1) ne $oligo)  ?  
          ( regres($oligo, 1), regres(complement($oligo, 1), 1) )  :  
          ( regres($oligo, 1) );
    $OL_DATA{REGEX}->{$oligo} = \@arr;            
  }
  #close IN;
  return \%OL_DATA;
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

=head2 make_oligos()
  takes a nucleotide sequence from a Chunk object and breaks it into oligos. A 
  hash reference provides all of the options, like target subchunk length, 
  oligo number, oligo length, etc. returns all oligos on one strand!!!!!!!!!!
  in: Chunk (struct), Options (hash reference)
  out: hashref

=head2 compare_sequences()
  takes two nucleotide sequences that are assumed to be perfectly aligned and 
  roughly equivalent and returns similarity metrics. should be twweeaakkeedd
  in: 2x nucleotide sequence (string)
  out: similarity metrics (hash)

=head2 count()
  takes a nucleotide sequence and returns a base count.  Looks for total length,
  purines, pyrimidines, and degenerate bases. If degenerate bases are present 
  assumes their substitution for non degenerate bases is totally random for 
  percentage estimation. 
  in: nucleotide sequence (string), 
  out: base count (hash)

=head2 ntherm()
  takes a nucleotide sequence and returns entropy, enthalpy, and free energy.   
  in: nucleotide sequence (string) 
  out: (array of integers) entropy enthalpy free energy 

=head2 compareseqs()
  takes nucleotide sequences and returns 1 if either could be said to be a
  perfect or degenerate copy of the other
  in: 2x nucleotide sequence (string)
  out: 1 OR 0

=head2 regres()
  takes a  sequence that may be degenerate and returns a string that is prepped
  for use in a regular expression.
  in: sequence (string), 
      switch for aa or nt sequence (1 or null)
  out: regexp string (string)

=head2 complement()
  takes a nucleotide sequence and returns its complement or reverse complement.
  in: nucleotide sequence (string), 
      switch for reverse complement (1 or null)
  out: nucleotide sequence (string)

=head2 melt()
  takes a nucleotide sequence and returns a melting temperature.  Has four 
  different formulas: 1 simple, 2 baldwin, 3 primer3, or 4 nntherm
  in: nucleotide sequence (string), 
      formula number, 
      salt concentration (string, opt, def =.05), 
      oligo concentration (string, opt, def = .0000001)
  out: temperature (string)

=head2 cleanup()
  takes a sequence and attempts to remove extraneous information.
  in: nucleotide sequence (string),
      switch for sequence type (0 strict nt, 1 degenerate nt, or 2 aa)
  out: nucleotide sequence  (string)

=head2 oligocruncher()
  takes a nucleotide sequence from a Chunk object and breaks it into oligos.  
  A hash reference provides all of the options, like target subchunk length, 
  oligo number, oligo length, etc   
  in: Chunk (struct)
      Options (hash reference)
  out: nothing (modifies Chunk (struct))

=head2 orf_finder()

=head2 define_oligos()

=head2 fasta_parser()

=head2 cons_seq()

=head2 print_alignment()
  $swit = 1 for html, 0 for text

=head2 fasta_writer()

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
