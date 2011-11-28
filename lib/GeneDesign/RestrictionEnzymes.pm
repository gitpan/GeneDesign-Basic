package GeneDesign::RestrictionEnzymes;
require Exporter;

use GeneDesign::Basic qw($ambnt);
use GeneDesign::SufTree;
use Perl6::Slurp;

use strict;
use warnings;

use vars qw($VERSION @ISA @EXPORT_OK %EXPORT_TAGS);
$VERSION = 3.00;

@ISA = qw(Exporter);
@EXPORT_OK = qw(
  define_sites
  overhang
  define_site_status
  siteseeker
  filter_sites
  mutexclu
  first_base
  report_RE
  build_suffix_tree
  search_suffix_tree
  search_DNA
  overhangP
  overhangA
  $IIA $IIA2 $IIA3 $IIP $IIP2 $IIPreg $IIAreg $IIBreg $treehit
);

%EXPORT_TAGS =  (all => [qw(define_sites overhang define_site_status siteseeker
  filter_sites mutexclu first_base report_RE build_suffix_tree
  search_suffix_tree search_DNA overhangP overhangA $IIA $IIA2 $IIA3 $IIP $IIP2
  $IIPreg $IIAreg $IIBreg $treehit)]);

our $IIA    = qr/ \( (    \d+) \/ (   \d+) \)  /x;
our $IIA2    = qr/ \( (    \d+) \/ (\- \d+) \)  /x;
our $IIA3    = qr/ \( (\-  \d+) \/ (\- \d+) \)  /x;
our $IIP    = qr/    ([A-Z]+) \^ [A-Z]*      /x;
our $IIP2    = qr/ \^ ([A-Z]+)          /x;
our $IIPreg  = qr/   ([A-Z]*)   \^ ([A-Z]*)      /x;
our $IIAreg  = qr/\A \w+ \(([\-]*\d+) \/ ([\-]*\d+)\)\Z  /x;
our $IIBreg  = qr/\A\(([\-]*\d+) \/ ([\-]*\d+)\) \w+ \(([\-]*\d+) \/ ([\-]*\d+)\)\Z  /x;
our $treehit  = qr/(\d+): ([A-Z0-9a-z]+) ([A-Z]+)/;

################################################################################
######################### Restriction Enzyme Functions #########################
################################################################################
sub define_sites
{
  my ($file) = @_;
  my @data = split(/\n/, slurp($file));
  my $swit = $data[0] =~ /^name/ ? 1  :  0;
  shift @data if ($swit == 1);
  my %RE_DATA;
  ($RE_DATA{REGEX},  $RE_DATA{CLEAN},  $RE_DATA{TABLE},  $RE_DATA{SCORE},  $RE_DATA{TYPE}, $RE_DATA{CLASS}, 
   $RE_DATA{DAM},    $RE_DATA{DCM},    $RE_DATA{CPG},    $RE_DATA{VEND} ,  $RE_DATA{STAR},  $RE_DATA{TEMP},
   $RE_DATA{BUF1},  $RE_DATA{BUF2},    $RE_DATA{BUF3} ,  $RE_DATA{BUF4} ,  $RE_DATA{BUFU}, $RE_DATA{INACT}, 
   $RE_DATA{EXCLUDE}    ) = 
  (  {}, {}, {}, {},   {}, {}, {},      {}, {}, {},      {}, {}, {},      {},  {}, {},    {}, {}, {});
  foreach my $line (@data)
  {
    my ($name, $site, $temp, $inact, $buf1, $buf2, $buf3, $buf4, $bufu, $dam, $dcm, $cpg, $score, $star, $vendor) = split("\t", $line);
    $RE_DATA{TABLE}->{$name}  = $site;    #annotated cleavage site
    $site =~ s/\W*\d*//g;
    $RE_DATA{CLEAN}->{$name}  = $site;    #recognition site
    $RE_DATA{TEMP} ->{$name}  = $temp;    #optimal temperature 
    $RE_DATA{INACT}->{$name}  = $inact;    #inactivation formula 
    $RE_DATA{BUF1} ->{$name}  = $buf1; 
    $RE_DATA{BUF2} ->{$name}  = $buf2; 
    $RE_DATA{BUF3} ->{$name}  = $buf3; 
    $RE_DATA{BUF4} ->{$name}  = $buf4; 
    $RE_DATA{BUFU} ->{$name}  = $bufu;
    $RE_DATA{DAM}  ->{$name}  = $dam; 
    $RE_DATA{DCM}  ->{$name}  = $dcm; 
    $RE_DATA{CPG}  ->{$name}  = $cpg;  
    $RE_DATA{SCORE}->{$name}  = $score;    #price of the enzyme
    $RE_DATA{STAR} ->{$name}  = $star;    #star activity?
    $RE_DATA{VEND} ->{$name}  = $vendor;    #vendors provided by
    my $sitelen = length($site);
    my ($lef, $rig) = ("", "");
    if ($RE_DATA{TABLE}->{$name} =~ $IIPreg)
    {
      $lef = length($1); 
      $rig = length($2);
      $RE_DATA{CLASS}->{$name} = "IIP";
    }
    elsif ($RE_DATA{TABLE}->{$name} =~ $IIBreg)
    {
      $lef = int($1);
      $rig = int($2);
      $RE_DATA{CLASS}->{$name} = "IIB";
    }
    elsif ($RE_DATA{TABLE}->{$name} =~ $IIAreg)
    {
      $lef = int($1); 
      $rig = int($2);
      $RE_DATA{CLASS}->{$name} = "IIA";
    }
  #stickiness
    $RE_DATA{TYPE}->{$name} .= ($lef < $rig)  ?  "5'"  :  ($lef > $rig)  ?  "3'"  :  ($lef == $rig)  ?  "b"  :  "?";
    $RE_DATA{TYPE}->{$name} .= "1" if (abs($lef - $rig) == 1);
  #regular expression array  
    my $arr = ( complement($site, 1) ne $site )            ?  
          [ regres($site, 1), regres(complement($site, 1), 1) ]  :  
          [ regres($site, 1) ];
    $RE_DATA{REGEX}->{$name} = $arr;            
  }
  foreach my $re (sort keys %{$RE_DATA{CLEAN}})
  {
    my %excl;
    foreach my $ar (sort grep {$_ ne $re} keys %{$RE_DATA{CLEAN}})
    {
      foreach my $arreg (@{$RE_DATA{REGEX}->{$ar}})
      {
        $excl{$ar}++ if ($RE_DATA{CLEAN}->{$re} =~ $arreg)
      }
      foreach my $rereg (@{$RE_DATA{REGEX}->{$re}})
      {
        $excl{$ar}++ if ($RE_DATA{CLEAN}->{$ar} =~ $rereg)
      }
    }
    my @skips = keys %excl;
    $RE_DATA{EXCLUDE}->{$re} = \@skips;
  }
  return \%RE_DATA;
}

sub report_RE
{
  my ($name, $RE_DATA) = @_;
  my $string = "$name: ", ;
  $string .= $$RE_DATA{TABLE}->{$name} . " (" . $$RE_DATA{CLEAN}->{$name} . 
      "), type = " . $$RE_DATA{TYPE}->{$name} ;
  return $string;
}

sub define_site_status
{
  my($seq, $SITE_REGEX) = @_;
  my $SITE_STATUS = {};
  foreach my $pattern (keys %$SITE_REGEX)
  { 
    my $count = 0;
    foreach my $sit (@{$$SITE_REGEX{$pattern}})
    {    
      $count++ while ($seq =~ /(?=$sit)/ig);
    }
    $$SITE_STATUS{$pattern} = $count;
  }
  return $SITE_STATUS;  
}

sub siteseeker
{
  my ($seq, $pass, $SITE_REGEX_arrref) = @_;
  my $total = {};
  foreach my $sit (@$SITE_REGEX_arrref)
  {
    while ($seq =~ /(?=$sit)/ig)
    {
      $$total{pos $seq} = substr($seq, pos $seq, length($pass));
    }
  }
  return $total;
}

sub build_suffix_tree
{
  my ($list_ref, $RE_DATA, $CODON_TABLE, $xlationref) = @_;
  my ($fehsh, $rehsh, $f_hsh, $r_hsh) = ({}, {}, {}, {});
  my ($ftree, $rtree) = (new_aa SufTree(), new_aa SufTree());
  foreach my $enzyme (@$list_ref)
  {
    my $site = $$RE_DATA{CLEAN}->{$enzyme};
    #    print "\t\t adding $site to tree\n";
    my $lagcheck = $site;
    $lagcheck = substr($lagcheck, 0, length($lagcheck)-1) while (substr($lagcheck, -1) eq "N");
    foreach my $peptide (amb_translation($site, $CODON_TABLE, $xlationref))
    {
     # $$fehsh{$peptide} = [] if (! exists $$fehsh{$peptide});
      push @{$$fehsh{$peptide}}, $enzyme;
    }
    my $etis = complement($site, 1);
    if ($etis ne $site && $lagcheck eq $site)
    {
      #      print "\t\t\t adding $etis to tree\n";
      foreach my $peptide (amb_translation($etis, $CODON_TABLE, $xlationref))
      {
        #  $$rehsh{$peptide} = [] if (! exists $$rehsh{$peptide});
        push @{$$rehsh{$peptide}}, $enzyme;      
      }
    }
  }
  $ftree->add_aa_paths($fehsh);
  $rtree->add_aa_paths($rehsh);
  return [$ftree, $rtree];
}

sub search_suffix_tree 
{
  my ($treelist, $exonlist, $RE_DATA, $CODON_TABLE, $xlationref) = @_;
  my $results;
  foreach my $exon (grep {$_} @$exonlist)
  {
    my $nucseq = $exon->seq;
       $nucseq = $nucseq->seq if ref $nucseq;
    my $exonname = $exon->display_name;
    my $aaseq = translate($nucseq, 1 + $exon->phase, $CODON_TABLE);
    #       print "Processing exon $exon ", $exon->display_name, " (", $exon->start, "..", $exon->end, ")\n";
 #   print "\t\t$nucseq<br>\n$aaseq\n\n";
    my $flip = 0;
    my $report = qr/\?/;
    my $reportenz = qr/\?/;
    foreach my $tree (@$treelist)
    {
      # print "\tSearching $flip tree...\n";
      $flip++;
      foreach my $rog ($tree->find_aa_paths($aaseq))
      {
        my $enz       = $$rog[0];
        my $recsite   = $$RE_DATA{CLEAN}->{$enz};
           $recsite   = complement($recsite, 1) if ($flip > 1);
        my $cutsite   = $$RE_DATA{TABLE}->{$enz};
        my $enzclass  = $$RE_DATA{CLASS}->{$enz};
        my $nucstart  = $$rog[1] * 3;
        my $sitestart = $exon->start + $nucstart - 1;
        my $sitestop  = $sitestart + (length($recsite)) - 1;
        my $peptide   = $$rog[2];
        my $presence  = "p";
        my $ohang     = {};
        my ($mattersbit, $ohangstart, $ohangend, $fabric, $offset, $situ, $pproof, $sitelen) = ("", 0, 0, "", 0, "", "",0);
        
        #                print "\tgot $enz @ $sitestart, $peptide\n";
      ##Figure Possible overhangs
        if ($$RE_DATA{TYPE}->{$enz} =~ /b/ || $$RE_DATA{CLASS}->{$enz} eq "IIB")
        {
          $situ = substr($nucseq, $nucstart, length($peptide)*3);
          $pproof = translate($situ, 1, $CODON_TABLE);
          $ohangstart = undef;
        }
        elsif ($$RE_DATA{TABLE}->{$enz} =~ $IIPreg)
        {
          my ($lef, $rig) = (length($1), length($2));
          ($lef, $rig) = ($rig, $lef) if ($rig < $lef);
          $ohangstart = length($recsite) - $rig + 1;
          $ohangend = length($recsite) - $lef;
          $situ = substr($nucseq, $nucstart, length($peptide)*3);
          ($fabric, $offset) = pattern_aligner($situ, $recsite, $peptide, $CODON_TABLE, 1, $xlationref);
          $mattersbit = substr($situ, $ohangstart + $offset -1, $ohangend - $ohangstart + 1);
          $pproof = translate($situ, 1, $CODON_TABLE);
          $sitelen = length($recsite);
        }
        elsif ($$RE_DATA{TABLE}->{$enz} =~ $IIAreg)
        {
          my ($lef, $rig) = ($1, $2);
          ($rig, $lef) = ($lef, $rig) if ($rig < $lef);
          $sitelen = $rig >= 0 ? length($recsite) + $rig  : length($recsite);
          my $nuclen = length($peptide)*3;
          $nuclen++ while($nuclen % 3 != 0);    
          $situ = substr($nucseq, $nucstart, $nuclen);
          ($fabric, $offset) = pattern_aligner($situ, $recsite, $peptide, $CODON_TABLE, 1, $xlationref);
          my $add;
          if ($flip == 1)
          {
            $ohangstart = length($recsite) + $lef + 1;
            unless ($rig <= 0)
            {
              $add = $rig - (length($fabric) - ($offset + length($recsite)));
              $add ++ while ($add % 3 != 0);
              $situ .= substr($nucseq, $nucstart + $nuclen, $add);
              $fabric .= "N" x $add;
            }
            
          }
          else
          {
            unless ($rig <= 0)
            {
              $add =  $rig - $offset;
              $add ++ while ($add % 3 != 0);
              $situ = substr($nucseq, $nucstart - $add, $add) . $situ;
              $fabric = "N" x $add . $fabric;
              $nucstart = $nucstart - $add;
              $ohangstart = $add - $rig + 1;
            }
            else
            {
              $ohangstart = $offset + abs($rig) + 1;
            }
          }
          $mattersbit = substr($nucseq, $nucstart + $ohangstart+1, $rig-$lef);
          $pproof = translate($situ, 1, $CODON_TABLE);
        }
        else
        {
          print "I don't recognize this type of enzyme: ";
          print "\t\t$enz $cutsite, $exonname\n\n";  
        }
        unless ($$RE_DATA{TYPE}->{$enz} =~ /b/ || $$RE_DATA{CLASS}->{$enz} eq "IIB")
        {
          if ($fabric eq "0")
          {
            print "oh no bad fabric, $enz, $fabric, $sitestart, $peptide, $exon\n";
            next; 
          }
          my $matstart = $ohangstart + $offset - 1;
             $matstart-- while($matstart % 3 != 0);
          my $matend = $ohangstart + $offset + length($mattersbit) - 1;
             $matend++ while($matend % 3 != 2);
          my $matlen = $matend - $matstart + 1;
          my $peproof = substr($pproof, ($matstart/3), $matlen/3);
          my $what = substr($fabric, $matstart, $matlen);
          #      print "\n\n$flip $enz $cutsite, $peptide, fab: $fabric, situ: $situ, pproof: $pproof, ohangstart: $ohangstart, ohangend: $ohangend, matters: $mattersbit, offset: $offset, matstart: $matstart, matlen: $matlen, peproof: $peproof, what: $what\n";
          #     print "\t\t $exonname: $nucstart of ", $exon->start, "..", $exon->end, "\n" ;#if ($matlen + $matstart >= $fabric);      
          foreach my $swapseq (amb_transcription($what, $CODON_TABLE, $peproof))
          {
            substr($fabric, $matstart, $matlen) = $swapseq;
            my $tohang = substr($fabric, $ohangstart +  $offset - 1, length($mattersbit));
            #      print "\t\t\tsubbing $swapseq into $fabric\tconsidering $tohang\n";
            $$ohang{$tohang}++ if ($tohang);#if ($tohang ne complement($tohang, 1));
          } 
          my @ohangers = keys %$ohang;
          #    print "\t\tpicked overhangs @ohangers\n";
        }
        
      ##Determine Presence
        $presence = "e" if ($situ =~ $$RE_DATA{REGEX}->{$enz}->[$flip - 1]);
        #        print "\t $enz $situ \n";
        $sitestart = $exon->start + $nucstart if ($exon->strand  != -1);
        $sitestart = $exon->end  - $nucstart + 1 - length($recsite) if ($exon->strand == -1);
        if ($sitestart + $sitelen - 1 > $exon->end)
        {
          $ohang = {$mattersbit => 1};
        }
        #    print "got $enz @ $sitestart, $presence, $peptide, $fabric\n" if ($presence eq "e");
        unless ($$RE_DATA{CLASS}->{$enz} eq "IIB" && $presence eq "p")
        {
          push @$results, [$enz, $sitestart, $presence, $ohang, $exon, $pproof, $ohangstart + $offset -1, $flip];
        }
      }
    }
  }
  return $results;
}

sub search_DNA
{
  my ($featurelist, $chrseq, $RE_DATA, $CODON_TABLE) = @_;
  my $results;
  foreach my $feature (@$featurelist)
  {
    #   print "Processing region ", $feature->display_id, " (", $feature->start, "..", $feature->stop, ")\n";
    my $start = $feature->start;
    my $nucseq = $feature->seq->seq;
    my $SITESTATUS = define_site_status($nucseq, $$RE_DATA{REGEX});
    foreach my $enz ( grep {$$SITESTATUS{$_} >= 1} keys %$SITESTATUS)
    {
      my $recsite = $$RE_DATA{CLEAN}->{$enz};
      my $regsite = $$RE_DATA{REGEX}->{$enz};
      my $tabsite = $$RE_DATA{TABLE}->{$enz};
      my $enzlocs = siteseeker($nucseq, $recsite, $regsite);
      foreach my $enzpos (keys %$enzlocs)
      {
        my $siteseq = $$enzlocs{$enzpos};
        my $strand = ($siteseq =~ $$RE_DATA{REGEX}->{$enz}->[0]) ? 1  : -1;
        my $sitestart = $enzpos + $start+1;
        my $presence = "i";
        my ($ohangoffset, $ohangseq) = (0, "");
        if ($$RE_DATA{CLASS}->{$enz} eq "IIP")
        {
          ($ohangoffset, $ohangseq) = overhangP($siteseq, $enz, $RE_DATA);
        }
        elsif ($$RE_DATA{TABLE}->{$enz} =~ $IIAreg)
        {
          my ($lef, $rig) = ($1, $2);
          ($lef, $rig) = ($rig, $lef) if ($rig < $lef);
          my $newseq;
          if ($strand == 1)
          {
            $newseq = substr($chrseq, $sitestart-2, length($recsite) + $rig + 5);
          }
          else
          {
            $newseq = substr($chrseq, $sitestart - ($rig+7), length($recsite) + $rig+5);
          }
          ($ohangoffset, $ohangseq) = overhangA($siteseq, $newseq, $enz, $RE_DATA, $strand);
        }
        elsif ($$RE_DATA{TABLE}->{$enz} =~ $IIBreg)
        {
          my ($rlef, $rrig) = ($3, $4);
          ($rlef, $rrig) = ($rrig, $rlef) if ($rrig < $rlef);
          ($ohangoffset, $ohangseq) = ($rlef, "NULL");
        }
        my $ohang = $ohangseq ?  {$ohangseq => 1} : {};
        push @$results, [$enz, $sitestart, $presence, $ohang, $feature, undef, $ohangoffset, $strand];
      }
    }    
  }
  return $results;
}

sub first_base
{
  my ($relev, $index, $orient) = @_;
  my $term = (($index -($relev %3))%3);
  if ($orient eq "+")
  {
    return $index - $term;
  }
  else
  {
    return $term != 0  ?  $index + (3-$term)  :  $index;
  }
}

sub overhang
{
  my ($dna, $pos, $grabbedseq, $table, $clean, $swit) = @_;
  my ($lef, $rig, $mattersbit, $cutoff) = (0, 0, undef, 0);
  my $orient = $grabbedseq =~ regres($clean)  ? "+"  :  "-";
  $swit = 0 if (!$swit);
  if ($table =~ $IIA || $table =~ $IIA2 || $table =~ $IIA3)
  {
    my ($lef, $rig) = ($1, $2);
    ($lef, $rig) =  (abs($lef), abs($rig)) if ($rig < 0 && $lef < 0);
    ($rig, $lef) = ($lef, $rig) if ($rig < $lef);
    $mattersbit = substr($dna, $pos + length($clean) + $lef, $rig - $lef) if ($orient eq "+");
    $mattersbit = substr($dna, $pos - $rig, $rig - $lef) if ($orient eq "-");
#          print "\t\t mattersbit: $mattersbit\n";
    $lef = length($grabbedseq) + $lef if ($lef < 0);
    return $swit != 1  ?  $mattersbit  :  [$mattersbit, $lef];
  }
  elsif ($table =~ $IIP || $table =~ $IIP2)
  {
    $cutoff = length($1);
    if ($cutoff > (.5 * length($clean))) { $cutoff = length($clean) - $cutoff;} 
    $mattersbit = substr($grabbedseq, $cutoff, length($grabbedseq)-(2*$cutoff));
    return $swit != 1  ?  $mattersbit  :  [$mattersbit, $cutoff];
  }
}

sub overhangP
{
  my ($dna, $enz, $RE_DATA) = @_;
  my ($ohangstart, $mattersbit) = (0, "");
  if ($$RE_DATA{TABLE}->{$enz} =~ $IIPreg)
  {
    my ($lef, $rig) = (length($1), length($2));
    ($lef, $rig) = ($rig, $lef) if ($rig < $lef);
    $ohangstart = $lef + 1;
    $mattersbit = substr($dna, $ohangstart-1, $rig-$lef);
    #    print "$enz, $dna, $context, ", $$RE_DATA{TABLE}->{$enz}, ", $mattersbit $ohangstart\n";
  }
  return ($ohangstart, $mattersbit);
}

sub overhangA
{
  my ($dna, $context, $enz, $RE_DATA, $strand) = @_;
  my ($ohangstart, $mattersbit) = (0, "");
  if ($$RE_DATA{TABLE}->{$enz} =~ $IIAreg)
  {
    my ($lef, $rig) = ($1, $2);
    ($lef, $rig) = ($rig, $lef) if ($rig < $lef);
    if ($strand == 1)
    {
      $ohangstart = length($dna) + $lef + 1;
    }
    else
    {
      $ohangstart = length($context) - length($dna) - $rig + 1;
    }
    $mattersbit = substr($context, $ohangstart-1, $rig-$lef);
    $ohangstart = $strand == 1  ? length($dna) + $lef :   0 - ($rig);
    #    print "$enz, $dna, $context, ", $$RE_DATA{TABLE}->{$enz}, ", $mattersbit $ohangstart ($strand)\n";
  }
  return ($ohangstart, $mattersbit);
}

sub filter_sites
{
  my ($pa_ref, $RE_DATA, $list_ref) = @_;
  my %pa = %$pa_ref;
  my @cutters = $list_ref  ?  @{$list_ref} : keys %{$$RE_DATA{CLEAN}};
  
  if ( $pa{check_stickiness} )
  {
    my %temp;
    my $regres = "[" . join("", split(" ", $pa{stickiness})) . "]";
    @cutters = grep {$$RE_DATA{TYPE}->{$_} =~ $regres} @cutters;
  }
  if ( $pa{check_cleavage_site} )
  {
    @cutters = grep { ($pa{cleavage_site} =~ /P/ && ($$RE_DATA{TABLE}->{$_} =~ $IIP || $$RE_DATA{TABLE}->{$_} =~ $IIP2))
            ||($pa{cleavage_site} =~ /A/ && ($$RE_DATA{TABLE}->{$_} =~ $IIA || $$RE_DATA{TABLE}->{$_} =~ $IIA2 || $$RE_DATA{TABLE}->{$_} =~ $IIA3)) }  @cutters;
  }
  if ( $pa{check_overhang} )
  {
    my %temp = map {$_ => 0} grep {$$RE_DATA{TYPE}->{$_} !~ /[1b]/} @cutters;
    foreach (keys %temp)
    {
      if ($$RE_DATA{TABLE}->{$_} =~ $IIP || $$RE_DATA{TABLE}->{$_} =~ $IIP2)
      {
        my $cutoff = length($1);
        my $clean = $$RE_DATA{CLEAN}->{$_};
        $cutoff = length($clean) - $cutoff if ($cutoff > (.5 * length($clean))); 
        my $mattersbit = substr($clean, $cutoff, length($clean)-(2*$cutoff));
        if ($mattersbit =~ $ambnt && length($mattersbit) % 2 == 0)
        {
          $temp{$_}++ if ($pa{overhang} =~ 'A');
        }
        else
        {
          $temp{$_}++ if (($mattersbit eq complement($mattersbit, 1) && length($mattersbit) % 2 == 0) && $pa{overhang} =~ 'P');
          $temp{$_}++ if (($mattersbit ne complement($mattersbit, 1) || length($mattersbit) % 2 == 1) && $pa{overhang} =~ 'N');
        }
      }
      $temp{$_}++ if (($$RE_DATA{TABLE}->{$_} =~ $IIA || $$RE_DATA{TABLE}->{$_} =~ $IIA2 || $$RE_DATA{TABLE}->{$_} =~ $IIA3) && $pa{overhang} =~ 'A');
    }
    @cutters = grep {$temp{$_} != 0} keys %temp;
  }
  if ( $pa{check_ambiguity} )
  {
    @cutters = grep {  ($pa{ambiguity} =~ /1/ && $$RE_DATA{CLEAN}->{$_} !~ /N/)
            || ($pa{ambiguity} =~ /2/ && $$RE_DATA{CLEAN}->{$_} !~ $ambnt) } @cutters;
  }
  if ( $pa{check_buffers})
  {
    if ($pa{buffer_bool} eq "OR")
    {
      my %temp;
      foreach my $b (split(" ", $pa{buffers}))
      {
        $temp{$_}++ foreach( grep { $$RE_DATA{"BUF".$b}->{$_} >= $pa{buffer_activity} } @cutters);
      }
      @cutters = keys %temp;
    }
    else
    {
      foreach my $b (split(" ", $pa{buffers}))
      {
        @cutters = grep { $$RE_DATA{"BUF".$b}->{$_} >= $pa{buffer_activity} } @cutters;
      }
    }
  }
  if ( $pa{check_heat})
  {
    my %temp;
    foreach my $h (split(" ", $pa{heat}))
    {
      $temp{$_}++ foreach( grep {$$RE_DATA{INACT}->{$_} =~ /\@$h/} @cutters);
    }
    @cutters = keys %temp;
  }
  if ( $pa{check_temperature})
  {
    my %temp;
    foreach my $t (split(" ", $pa{temperature}))
    {
      $temp{$_}++ foreach( grep {$$RE_DATA{TEMP}->{$_} == $t} @cutters);
    }
    @cutters = keys %temp;
  }
  if ( $pa{check_star} )
  {
    @cutters = grep {  ($pa{check_star} =~ /1/ && $$RE_DATA{STAR}->{$_} eq "y")
            || ($pa{check_star} =~ /2/ && !$$RE_DATA{STAR}->{$_}) } @cutters;
  }
  if ( $pa{check_price} )
  {
    @cutters = grep { $$RE_DATA{SCORE}->{$_} >= $pa{low_price} && $$RE_DATA{SCORE}->{$_} <= $pa{high_price} } @cutters;
  }
  if ( $pa{check_site_length} )
  {
    my %lenged;
    foreach my $len ( split " ", $pa{site_length} )
    {
      if ($len ne 'b')
      {
        $lenged{$_}++ foreach ( grep {length($$RE_DATA{CLEAN}->{$_}) == $len} @cutters);
      }
      else
      {
        $lenged{$_}++ foreach ( grep {length($$RE_DATA{CLEAN}->{$_}) > 8} @cutters);
      }
    }
    @cutters = keys %lenged;
  }
  if ( $pa{check_meth_status} )
  {
    if ($pa{check_meth_status} == 1)
    {
      my %methed;
      $methed{$_}++ foreach ( grep { ! $$RE_DATA{DAM}->{$_} &&  ! $$RE_DATA{CPG}->{$_}  && ! $$RE_DATA{DCM}->{$_} } @cutters);
      @cutters = keys %methed;
    }
    elsif ($pa{check_meth_status} == 2)
    {
      if ($pa{meth_bool} eq "OR")
      {
        my %temp;
        foreach my $methstat (split " ", $pa{meth_status})
        {
          my ($type, $way) = ($1, uc $2) if ($methstat =~ /(\w)\.(\w+)/);
          $temp{$_}++ foreach( grep { $$RE_DATA{$way}->{$_} eq $type || ($type eq "f" && ! $$RE_DATA{$way}->{$_}) } @cutters);
        }
        @cutters = keys %temp;
      }
      else
      {
        foreach my $methstat (split " ", $pa{meth_status})
        {
          my ($type, $way) = ($1, uc $2) if ($methstat =~ /(\w)\.(\w+)/);
          @cutters = grep { $$RE_DATA{$way}->{$_} eq $type || ($type eq "f" && ! $$RE_DATA{$way}->{$_}) } @cutters;
        }
      }
    }
  }
  foreach my $forbid (split " ", $pa{disallowed_seq})
  {
    @cutters = grep { $$RE_DATA{CLEAN}->{$_}  !~ regres($forbid) } @cutters;
  }
  my %forced;
  foreach my $force (split " ", $pa{required_seq})
  {
    $forced{$_}++ foreach (grep { $$RE_DATA{CLEAN}->{$_}  =~ regres($force) } @cutters);
    @cutters = keys %forced;
  }

  return @cutters;
}

sub mutexclu
{
  my ($used_ref, $allsites_ref) = @_;
  #don't want second degree exclusion - only exclude if it's not a -2 site
  foreach my $c ( grep {$$used_ref{$_} != -2} keys %{$used_ref})
  {
    my $te = regres($$allsites_ref{$c});
    foreach my $d ( grep {$c ne $_} keys %{$allsites_ref})
    {
      my $ue = regres($$allsites_ref{$d});
      $$used_ref{$d} = -2  if ($$used_ref{$c} =~ $ue || $$allsites_ref{$d} =~ $te);
    }
  }
  return $used_ref;
}

1;

__END__

=head1 NAME

GeneDesign::RestictionEnzymes

=head1 VERSION

Version 3.00

=head1 DESCRIPTION

  GeneDesign functions for handling restriction enzyme recognition sites
  
=head1 Functions

=head2 define_sites()
  Generates a set of hash references containing information about the 
  restriction enzyme library. All have the enzyme name as a key. 
    REGEX : ref to array containing regex-ready representations of enzyme 
            recognition site
    CLEAN : the recognition site as a string
    TABLE : the recognition site with cleavage annotation
    SCORE : usually the price
    TYPE  : whether the enzyme leaves 5', 3', or blunt ends, and whether or not 
            the enzyme leaves a 1bp overhang
    METHB : which methylation pattern blocks cleavage
    METHI : which methylation pattern inhibits cleavage
    STAR  : whether or not the enzyme has star activity
    TEMP  : optimal operating temperature for cleavage
    VEND  : which vendors supply the enzyme
  in: enzyme list as a file path
  out: hash of hash references

=head2 overhang()

=head2 define_site_status()
  Generates a hash describing a restriction count of a nucleotide sequence.
  in: nucleotide sequence as a string
      reference to a hash containing enzyme names as keys and regular
        expressions in a reference to an array as values
        (SITE_REGEX from define_sites)
  out: reference to a hash where the keys are enzyme names and the value is a 
        count of their occurence in the nucleotide sequence

=head2 siteseeker()
  Generates a hash describing the positions of a particular enzyme's recognition
  sites in a nucleotide sequence.
  in: nucleotide sequence as a string,
      an enzyme recognition site as a string, 
      reference to an array containing enzyme regular expressions
  out: reference to a hash where the keys are positions and the value is the 
      recognition site from that position as a string.

=head2 filter_sites()

=head2 mutexclu()

=head2 first_base()

=head2 report_RE()

=head2 build_suffix_tree()
  Builds two suffix trees from a list of restriction enzyme recognition sites
  in: array of enzyme names,
      reference enzyme hash
  out: array of suffix trees, one for forward orientaion, one for reverse

=head2 search_suffix_tree()
  Working on a Bio::Seq::Feature::Store DB

=head2 search_DNA()

=head2 overhangP()

=head2 overhangA()

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
