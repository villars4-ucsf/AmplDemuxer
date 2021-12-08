    cutadapt \
      --action=retain \
      --discard-untrimmed \
      -g fasta_files/v3_fwd_primers_1B_anchored.fasta \
      -G fasta_files/v3_rev_primers_1B_anchored.fasta \
      --pair-adapters \
      -o untrimmed_R1.fastq.gz \
      -p untrimmed_R2.fastq.gz \
      -e 0 \
      --no-indels \
      seq_files_test/test_S100_R1_001.fastq.gz \
      seq_files_test/test_S100_R2_001.fastq.gz \
      > cutadapt.log
