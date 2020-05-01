#!/bin/bash

#Treat expanding unset parameters as error 
set -u
#Exit with status 1 if a command throws an error
set -e

if [ $# -le 0 ]
then
  echo ""
  echo "Not enough arguments supplied. Execute as"
  echo "PATH_TO_RELATE/scripts/clues_Ag.sh"
  echo ""
  echo "-p,--position:  Position of SNP."
  echo "-n,--name:  Name of SNP." 

  exit 1;
fi


POSITIONAL=()
while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
    -p|--position)
      pos="$2"
      shift # past argument
      shift # past value
      ;;
    -n|--name)
      name="$2"
      shift # past argument
      shift # past value
      ;;
    *)    # unknown option
      POSITIONAL+=("$1") # save it in an array for later
      shift # past argument
      ;;
  esac
done


echo "Parameters passed to script:"
echo "pos         = $pos"
echo "pops        = $pops"
echo "name        = $name"


./home/sanj/apps/Relate/scripts/SampleBranchLengths/SampleBranchLengths.sh -i $pop -o "${pop}.${name}" -m 5.5e-9 --coal {params.coal} --num_samples {params.nsamples} --dist {params.dist} --first_bp $pos --last_bp $pos

