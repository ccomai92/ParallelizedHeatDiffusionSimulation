#!/usr/bin/env bash

# UWB MPI Cluster Setup Script
# Thomas Kercheval, 2019
#
# Configures SSH access to all machines in MPI cluster, then
#   sets user's password on each machine. The MPI cluster is
#   backed by a NFS server, so every student's home directory
#   is shared among machines. Changing `~/.ssh` on one machine
#   will propagate the change to all other servers.


if [ $# -ne 2 ]; then
    echo "ERROR: Invalid number of arguments"
    echo "Usage: ./setup_mpi_cluster.sh OLD_PASS NEW_PASS"
    exit
fi

# Number of client computers in cluster
NUM_CLIENTS=8

OLD_PASS="$1"
NEW_PASS="$2"
PASS_STRING="${OLD_PASS}\n${NEW_PASS}\n${NEW_PASS}"
SSH_PASS_STR="echo -e \"${PASS_STRING}\" | passwd"

# Create ~/.ssh file, if it does not exist
if [[ ! -d "$HOME/.ssh" ]]; then
  mkdir $HOME/.ssh
fi

chmod 700 $HOME/.ssh

# Trust all host fingerprints in cluster
for ((i=1; i <= $NUM_CLIENTS; i++)); do
  ssh-keyscan -H cssmpi$i.uwb.edu >> $HOME/.ssh/known_hosts;
done;

# Generate SSH key
ssh-keygen -N "" -f $HOME/.ssh/id_rsa 

# Add generated SSH key to trusted key list
cat $HOME/.ssh/id_rsa.pub >> $HOME/.ssh/authorized_keys
chmod 600 $HOME/.ssh/authorized_keys

# Set new password on each client
for ((i=1; i <= $NUM_CLIENTS; i++)); do
  ssh -i $HOME/.ssh/id_rsa $USER@cssmpi$i.uwb.edu ${SSH_PASS_STR};
done;
