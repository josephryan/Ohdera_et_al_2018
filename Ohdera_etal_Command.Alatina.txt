############## MaSuRCA assembly of Alatina

/usr/local/MaSuRCA-3.2.2/bin/masurca sr_config_example.txt --skip-checking
./assemble.sh

#orig assemble.sh ran for a while but then died cuz of too many openfiles
# had to up the amount of files available.

# http://posidev.com/blog/2009/06/04/set-ulimit-parameters-on-ubuntu/
# /etc/pam.d/pam_limits.conf  # for scientific linux vs ubuntu`

export LD_LIBRARY_PATH=/usr/local/MaSuRCA-3.2.2/lib
/usr/local/MaSuRCA-3.2.2/bin/masurca sr_config_example.txt --skip-checking
./assemble.sh

# failed again so then ran:
mv CA.mr.41.15.17.0.029/5-consensus/ ./99-FAILED_BITS/

export LD_LIBRARY_PATH=/usr/local/MaSuRCA-3.2.2/lib
/usr/local/MaSuRCA-3.2.2/bin/masurca sr_config_example.txt --skip-checking
./assemble.sh
