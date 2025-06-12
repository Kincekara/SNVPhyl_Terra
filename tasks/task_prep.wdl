version 1.0

task validate_inputs {
  input {
    Array[String] samplename
    Array[File] read1
    Array[File] read2
    String taxon
    File? reference    
    String? accession
  }

  command <<<
    set -e

    # Initialize INPUT.CHECK
    echo "PASS" > INPUT.CHECK

    # Convert WDL arrays to shell arrays    
    samples=(~{sep=' ' samplename})
    r1=(~{sep=' ' read1})
    r2=(~{sep=' ' read2})
    len=${#samples[@]}

    # Check if the lengths of read1 and read2 match the length of samplename
    if [ "${#r1[@]}" -eq "$len" ] && [ "${#r2[@]}" -eq "$len" ]; then
      for (( i=0; i<len; i++ )); do
        f1=$(basename "$r1[i]")
        f2=$(basename "$r2[i]")
        if [ "$f1" = "$f2" ]; then
          echo "Same read1 and read2 provided for $samples[i]"
          echo "FAIL" > INPUT.CHECK
        fi        
      done
      # Check INPUT.CHECK status before proceeding
      if [ "$(cat INPUT.CHECK)" = "FAIL" ]; then
        exit 1
      fi
      # Create samples.tsv
      for (( i=0; i<len; i++ )); do
        read1=$(echo "${r1[i]}" | sed 's/\/mnt\/disk\/cromwell_root/gs\:\//g')
        read2=$(echo "${r2[i]}" | sed 's/\/mnt\/disk\/cromwell_root/gs\:\//g')
        echo -e "${samples[i]}\t$read1\t$read2" >> "samples.tsv"
      done
    else
      echo "Missing reads file!"
      echo "FAIL" > INPUT.CHECK
      exit 1
    fi   

    # Check reference genome
    if [ -f "~{reference}" ]; then
      mv ~{reference} reference.fasta
    elif [ -n "~{accession}" ]; then
      datasets download genome accession ~{accession} --include genome --filename reference.zip
      unzip -j -p reference.zip *.fna > reference.fasta
    elif [ -n "~{taxon}" ]; then
      id=$(datasets summary genome taxon "~{taxon}" --reference | cut -d "," -f1 | cut -d ":" -f3 | sed 's/\"//g')
      datasets download genome accession "$id" --include genome --filename reference.zip
      unzip -j -p reference.zip *.fna > reference.fasta
    else
      echo "No reference information provided! Please give one of them: reference file, toxon or accession id"
      echo "FAIL" > INPUT.CHECK
      exit 1
    fi  
  >>>

  output {
    File ref_genome = "reference.fasta"
    File samplelist = "samples.tsv"
    String check = read_string("INPUT.CHECK")
  }

  runtime {
    docker: "staphb/ncbi-datasets:15.31.1"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}