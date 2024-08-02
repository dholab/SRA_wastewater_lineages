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

    # fetch the SRA accession
    echo $SRA_NUM
    python3 ./SRA_fetch.py --SRA=${SRA_NUM}

    # run sam_refiner
    python3 ./SAM_Refiner.py -r SARS2.gb \
    --wgs 1 --collect 0 --seq 1 --indel 0 --covar 0 --nt_call 1 --min_count 1 \
    --min_samp_abund 0 --ntabund 0 --ntcover 1 --mp 4 --chim_rm 0 --deconv 0  -S ${SRA_NUM}.SARS2.wg.sam

    # sort and convert the sam_refiner output to a CRAM
    samtools sort ${SRA_NUM}.SARS2.wg.sam \
    | samtools view -T SARS2.fasta -@8 -o ${SRA_NUM}.SARS2.wg.cram -

    # extract the unique variants
    python3 ./Variant_extractor.py ${SRA_NUM}

done
