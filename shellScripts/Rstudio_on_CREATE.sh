#!/bin/sh
#SBATCH --job-name=rstudio-create
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --mem=25G
#SBATCH --time=10:00:00

PATH=$R_HOME/bin:$PATH

module load rstudio_kcl/v2022.07.1_554-gcc-9.4.0-r4.1.1-python-3.8.12

export PASSWORD=$(openssl rand -base64 15)
export IP_ADD=$(hostname -I | awk '{print $1}')

# get unused socket per https://unix.stackexchange.com/a/132524
export PASSWORD=$(openssl rand -base64 15)
readonly IPADDRESS=$(hostname -I | tr ' ' '\n' | grep '10.211.4.')
readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
cat 1>&2 <<END
1. SSH tunnel from your workstation using the following command:

 

   ssh -NL 8787:${HOSTNAME}:${PORT} ${USER}@hpc.create.kcl.ac.uk

 

   and point your web browser to http://193.61.76.38:8787

 


2. Login to RStudio Server using the following credentials:

 

   user: ${USER}
   password: ${PASSWORD}

 

When done using the RStudio Server, terminate the job by:

 

1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:

 

   scancel -f ${SLURM_JOB_ID}

 

END

# Create custom database config
# $TMPDIR is where the database file and server data dir are stored - can amend to a different directory (e.g. /scratch/users/kXXXXXXXX)
DBCONF=$TMPDIR/database.conf
if [ ! -e $DBCONF ]
then
  printf "\nNOTE: creating $DBCONF database config file.\n\n"
  echo "directory=$TMPDIR/var-rstudio-server" > $DBCONF
fi

 

rserver --server-user ${USER} --www-port ${PORT} --secure-cookie-key-file $TMPDIR/secure-cookie-key --session-rpc-key-file=$TMPDIR/session-rpc-key  --server-data-dir $TMPDIR/data-rstudio-server --database-config-file=$DBCONF --auth-none=0 --auth-pam-helper-path=pam-env-helper

 

printf 'RStudio Server exited' 1>&2