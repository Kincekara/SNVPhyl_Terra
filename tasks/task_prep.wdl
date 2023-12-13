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
    # check sample inputs
    echo "PASS" > INPUT.CHECK    
    samples=(~{sep=' ' samplename})
    r1=(~{sep=' ' read1})
    r2=(~{sep=' ' read2})
    len=${#samples[@]}
    if [ ${#r1[@]} -eq $len ] && [ ${#r2[@]} -eq $len ]; then
      for (( i=0; i<$len; i++ )); do
        f1=$(basename $r1[i])
        f2=$(basename $r2[i])
        if [ $f1 == $f2]; then
          echo "Same read1 and read2 provided for $samples[i]"
          echo "FAIL" > INPUT.CHECK
        fi        
      done
      # create samples.tsv
      for (( i=0; i<$len; i++ )); do
        echo -e "${samples[i]}\t${r1[i]}\t${r2[i]}" >> "samples.tsv"
      done
    else
      echo "Missing reads file!"
      echo "FAIL" > INPUT.CHECK
    fi   

    # check reference genome
    if [ -f "~{reference}" ]; then
      mv ~{reference} reference.fasta
    elif [ -z "~{accession}" ]; then
      datasets download genome accession ~{accession} --include genome --filename reference.zip
      unzip -j -p reference.zip *.fna > reference.fasta
    elif [ -z "~{taxon}" ]; then
      id=$(datasets summary genome taxon "~{taxon}" --reference | cut -d "," -f1 | cut -d ":" -f3 | sed 's/\"//g')
      datasets download genome accession $id --include genome --filename reference.zip
      unzip -j -p reference.zip *.fna > reference.fasta
    else
      echo "No reference information provided! Please give one of them: reference file, toxon or accession id"
      echo "FAIL" > INPUT.CHECK
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