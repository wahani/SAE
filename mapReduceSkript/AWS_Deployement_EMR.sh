#!/bin/sh
# turn on logging and exit on error
set -e -x

# update packages
sudo apt-get update
sudo DEBIAN_FRONTEND=noninteractive apt-get upgrade --yes --force-yes

# altes R entfernen
sudo apt-get remove r-base r-cran-* --yes
sudo apt-get autoremove --yes

# add key
sudo apt-key adv --keyserver subkeys.pgp.net --recv-key 381BA480

# add cran-mirror to sources
sudo su -- -c 'echo "deb http://cran.r-project.org/bin/linux/debian squeeze-cran3/" >> /etc/apt/sources.list'

# pin cran-mirror to install recent version of R
sudo su -- -c 'echo "" >> /etc/apt/preferences'
sudo su -- -c 'echo "Package: *" >> /etc/apt/preferences'
sudo su -- -c 'echo "Pin: origin cran.r-project.org" >> /etc/apt/preferences'
sudo su -- -c 'echo "Pin-Priority: 1000" >> /etc/apt/preferences'

# install R and Git
sudo apt-get update
sudo apt-get install --yes r-recommended
sudo apt-get install --yes git
sudo apt-get clean

# optional project specific actions
echo -n "Starting additional bootstrap actions for project SAE... "

# pull the project
cd ~
git clone https://github.com/wahani/SAE.git
cd SAE
git pull origin mapReduce
git checkout mapReduce

# prepare R
cd ~
sudo su -- -c "R -e \"install.packages(c('testthat', 'Rcpp', 'RcppArmadillo'), repos='http://cran.rstudio.com/')\""
sudo R CMD INSTALL SAE/spatioTemporalData_1.1.5.tar.gz
sudo R CMD INSTALL --no-multiarch --with-keep.source SAE/package
echo "done."