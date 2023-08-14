tar -xzf R413.tar.gz
# (optional) if you have a set of packages (created in Part 1), untar them also
tar -xzf R_packages.tar.gz
export PATH=$PWD/R/bin:$PATH
export RHOME=$PWD/R
export R_LIBS=$PWD/packages
Rscript SBM.R 200
cat SBM.R
R
nano SBM.R
Rscript SBM.R 200
exit
