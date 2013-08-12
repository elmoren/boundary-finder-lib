use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'BoundaryFinder',
    VERSION_FROM => 'lib/BoundaryFinder.pm', # finds \$VERSION
    AUTHOR       => 'Nathan Elmore (nate@elmoren.com)',
    ABSTRACT     => 'Exon/Intron Boundary Finder',
    PREREQ_PM    => {
                     'File::Basename' => 0.01,
                     'File::Which' => 0.01,
                     'Test::Simple' => 0.44,
                    },

);
