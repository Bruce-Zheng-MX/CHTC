ls
tar -xzf R413.tar.gz
# (optional) if you have a set of packages (created in Part 1), untar them also
tar -xzf R_packages.tar.gz
export PATH=$PWD/R/bin:$PATH
export RHOME=$PWD/R
export R_LIBS=$PWD/packages
Rscript -e 'install.packages(c("wbs","not"), repos="https://cloud.r-project.org")'
Rscript real_data_application.R 100 1
cat real_data_application.R 
R
nano real_data_application.R 
ls
exit
