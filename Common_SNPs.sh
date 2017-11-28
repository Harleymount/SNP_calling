#!/bin/bash

#put reference in FA format with your read 1 and read 2  fastq in a directory that you provide with --directory 
#requirements novoalign, 2bwt-builder, SOAP, Picard, FreeBayes, 
for arg in "$@"; do
  shift
  case "$arg" in
    "--vcf_1") set -- "$@" "-r" ;;
    "--vcf_2") set -- "$@" "-x" ;;
	"--directory")   set -- "$@" "-z" ;;
    *)        set -- "$@" "$arg"
  esac
done


while getopts r:x:z: option
do
 case "${option}"
 in
 r) vcf1=${OPTARG};;
 x) vcf2=${OPTARG};;
 z) dir=$OPTARG;;
 esac
done

cd $dir;

Rscript common_SNP.R $vcf1 $vcf2;
