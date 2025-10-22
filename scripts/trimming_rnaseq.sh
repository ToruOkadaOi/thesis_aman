for sample in S55 S56 S57 S58 S59 S60; do
  if [[ "$sample" == "S59" ]]; then
    r1=fixed_S59_R1.fastq.gz
    r2=fixed_S59_R2.fastq.gz
  elif [[ "$sample" == "S60" ]]; then
    r1=fixed_S60_R1.fastq.gz
    r2=A006850256_186642_S60_L002_R2_001.fastq.gz
  else
    r1=$(ls A006850256_*_${sample}_L002_R1_001.fastq.gz)
    r2=$(ls A006850256_*_${sample}_L002_R2_001.fastq.gz)
  fi

  fastp -i "$r1" -I "$r2" \
        -o trimmed/${sample}_R1_trimmed.fastq.gz \
        -O trimmed/${sample}_R2_trimmed.fastq.gz \
        -j fastp_logs/${sample}_fastp.json \
        -h fastp_logs/${sample}_fastp.html \
        -w 8
done
