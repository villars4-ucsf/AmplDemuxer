version 1.0

workflow AmplDemuxer {
  input {
    String sample_name
    Array[File] read_samples
    Array[Array[File]] primer_sets
  }
  
  scatter(primer_set in primer_sets){
    call BinByAmplicon {input: read_samples = read_samples, primer_set = primer_set}
    call PrimerRemoval {input: binned_forward_seqs = BinByAmplicon.binned_forward_seqs, binned_reverse_seqs = BinByAmplicon.binned_reverse_seqs}
  }

  #call AmpliconDataGeneration {input: binned_forward_seqs_no_primer = PrimerRemoval.binned_forward_seqs_no_primer, binned_reverse_seqs_no_primer = PrimerRemoval.binned_reverse_seqs_no_primer}

  output {
    # Array[File] r2s = BinByAmplicon.r2
  }
}

task BinByAmplicon {
  input {
    Array[File] read_samples
    Array[File] primer_set
    }

  command {

    : '
      NOTE
        a. -j flag sets num cores
        b. primer seqs have caret in order to force primer search at begining of read
        c. primers need to be removed from reads bc the error
           profile of the primers skews the dada2 algorithm
    '

    # Save path to file containing path to primers
    primer_file=${write_lines(primer_set)}
    forward_primer=$(sed -n '1p' $primer_file)
    reverse_primer=$(sed -n '2p' $primer_file)

    # Save path to file containing path to sequences
    sequence_file=${write_lines(read_samples)}
    forward_reads=$(sed -n '1p' $sequence_file)
    reverse_reads=$(sed -n '2p' $sequence_file)

    # Removing Primers
    cutadapt \
      --action=retain \
      --discard-untrimmed \
      -g file:$forward_primer \
      -G file:$reverse_primer \
      --pair-adapters \
      -o untrimmed_R1.fastq.gz \
      -p untrimmed_R2.fastq.gz \
      -e 0 \
      --no-indels \
      $forward_reads \
      $reverse_reads > cutadapt.log

  }

   runtime {
     simg: "/wynton/home/rodriguez-barraquer/villars4_ucsf/library/singularity/cutadapt-2.5.sif"
   }

  output {
    File binned_forward_seqs = "binned_R1.fastq.gz"
    File binned_reverse_seqs = "binned_R2.fastq.gz"
    String log = "cutadapt.log"
  }
}

task PrimerRemoval {
  input {
    File binned_forward_seqs
    File binned_reverse_seqs
    }

  command {
    : '
      NOTE
        a. the "-w" flag sets thread count 
    '

    fastp \
      -i ${binned_forward_seqs} \
      -I ${binned_reverse_seqs} \
      -o binned_forward_seqs_no_primer.fastq.gz \
      -O binned_reverse_seqs_no_primer.fastq.gz \
      -h noprimers.html
  }

  # runtime {
  #   docker: "broadinstitute/my_image"
  # }

  output {
    File binned_forward_seqs_no_primer = "binned_forward_seqs_no_primer.fastq.gz"
    File binned_reverse_seqs_no_primer = "binned_reverse_seqs_no_primer.fastq.gz"
    File log = "noprimers.html"
  }
}

task AmpliconDataGeneration {
  input {
    File binned_forward_seqs_no_primer
    File binned_forward_seqs_no_primer
    }

  command {
    : '
      NOTE
        a. the "-w" flag sets thread count 
    '

    Rscript /Users/severianovillarruel/eppicenter/library/WDL/AmplDemuxer/dependant_scripts/test.R
  }

  # runtime {
  #   docker: "broadinstitute/my_image"
  # }

  output {

  }
}
