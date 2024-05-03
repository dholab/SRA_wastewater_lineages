sra_tracking_dir=/path/to/sra_tracking_dir
queued_filepath='${sra_tracking_dir}queued_list.txt'
for SRA_NUM in `cat ${sratxt}`;  do
  ./sra_cryptic_loop.sh
done