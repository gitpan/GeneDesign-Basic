use 5.006;
use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'GeneDesign::Basic',
    license             => 'bsd',
    dist_author         => q{Sarah Richardson <notadoctor@jhu.edu>},
    dist_version_from   => 'lib/GeneDesign/Basic.pm',
    build_requires => {
        'Test::More' => 0,
    },
    requires => {
        'perl' => 5.006,
        'List::Util' => 0,
        'Perl6::Slurp' => 0.03,
        'Text::Wrap' => 2009.0305,
    },
    configure_requires => {
      'Module::Build' => 0.38
    },
    add_to_cleanup      => [ 'GeneDesign-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
