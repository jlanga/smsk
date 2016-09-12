brew update

# To make the example work
brew install \
    python3 \
    graphviz

# Homebrew science packages to run the example
brew tap homebrew/science
brew install \
    homebrew/science/bwa \
    homebrew/science/samtools \
    homebrew/science/bcftools


# Examples of other languages: python3, R, ruby, perl
# brew install \
#     python3 \
#     homebrew/science/r \
#     ruby \
#     perl
