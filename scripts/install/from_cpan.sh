#!/usr/bin/env sh

# Install perl from brew
# brew install \
#    perl

# Download cpanminus (executable to install locally cpan modules, no sudo)
wget \
    --continue \
    --output-document bin/cpanm \
    https://raw.githubusercontent.com/miyagawa/cpanminus/master/cpanm
chmod +x bin/cpanm

# cpan
cpanm \
    Bio::Perl \
    Bit::Vector \
    DBD::SQLite \
    DBI \
    IO::All \
    IO::Prompt \
    Inline::C \
    Perl::Unsafe::Signals \
    forks \
    forks::shared
