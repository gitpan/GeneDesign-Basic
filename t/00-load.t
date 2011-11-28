#!perl -T

use Test::More tests => 6;

BEGIN {
    use_ok( 'GeneDesign::Basic' ) || print "Bail out!\n";
    use_ok( 'GeneDesign::Codons' ) || print "Bail out!\n";
    use_ok( 'GeneDesign::HTML' ) || print "Bail out!\n";
    use_ok( 'GeneDesign::Random' ) || print "Bail out!\n";
    use_ok( 'GeneDesign::RestrictionEnzymes' ) || print "Bail out!\n";
    use_ok( 'GeneDesign::SufTree' ) || print "Bail out!\n";
}

diag( "Testing GeneDesign::Basic $GeneDesign::Basic::VERSION, Perl $], $^X" );
