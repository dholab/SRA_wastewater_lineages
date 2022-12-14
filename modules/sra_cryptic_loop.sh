while getopts s: opt ; do
   case $opt in
      s) sratxt=$OPTARG ;;
      *) usage; exit 1;;
   esac
done

SRA_NUM=${sratxt%????}
echo ${SRA_NUM}
echo ${sratxt}
cat ${sratxt}
for SRA_NUM in `cat ${sratxt}`;  do
# rm -f ${sratxt}
  echo $SRA_NUM
  python3 ./SRA_fetch.py --SRA=${SRA_NUM}

  rm -rf ${SRA_NUM}
  rm -f *.merge.fq
  rm -f *.un1.fq
  rm -f *.un2.fq
  rm -f *_1.fastq
  rm -f *_2.fastq
  rm -f *.collapsed.fa

  python3 ./SAM_Refiner.py -r SARS2.gb \
  --wgs 1 --collect 0 --seq 1 --indel 0 --covar 0 --nt_call 1 --min_count 1 \
  --min_samp_abund 0 --ntabund 0 --ntcover 1 --mp 4 --chim_rm 0 --deconv 0  -S ${SRA_NUM}.SARS2.wg.sam

  samtools sort -o ${SRA_NUM}.SARS2.wg.srt.sam ${SRA_NUM}.SARS2.wg.sam
  samtools view -T SARS2.fasta -@8 ${SRA_NUM}.SARS2.wg.srt.sam -o ${SRA_NUM}.SARS2.wg.cram
  python3 ./Variant_extractor.py ${SRA_NUM}

  # mv dir_29618_AA_E484del.tsv ${SRA_NUM}_dir_29618_AA_E484del.tsv

  rm -f ${SRA_NUM}.SARS2.wg.sam
  rm -f ${SRA_NUM}.SARS2.wg.srt.sam
  gzip ${SRA_NUM}.SARS2.wg_unique_seqs.tsv
  rm -f ${SRA_NUM}.SARS2.wg_unique_seqs.tsv
  gzip ${SRA_NUM}.SARS2.wg_nt_calls.tsv
  rm -f ${SRA_NUM}.SARS2.wg_nt_calls.tsv
done

rm -f ./Variant_extractor.py
rm -f ./SRA_fetch.py
rm -f ./sra_cryptic.sh
rm -f ./SAM_Refiner.py
rm -f ./derep.py
rm -f ./SARS2.fasta
rm -f ./SARS2.gb
rm -f *.all.fq