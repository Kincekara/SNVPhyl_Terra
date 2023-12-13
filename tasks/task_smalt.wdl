version 1.0

task index {
  input {
    File genome
    String filename = basename(genome)    
  }

  command <<<
    smalt index -k 13 -s 6 ~{filename} ~{genome}
    cp ~{genome} ./~{filename}
    samtools faidx ~{filename}
  >>>

  output {
    File reference = "~{filename}"
    File sma = "~{filename}.sma"
    File smi = "~{filename}.smi"
    File fai = "~{filename}.fai"
  }

  runtime {
    docker: "staphb/smalt:0.7.6"
    memory: "2 GB"
    cpu: 1
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 1
  }
}

task map {
  input {
    String samplename
    File read1
    File read2    
    File genome        
    File fai
    File sma
    File smi
  }

  command <<<
    # map reads
    smalt map -f bam -n 4 -l pe -i 1000 -j 20 -r 1 -y 0.5 -o ~{samplename}.bam ~{genome} ~{read1} ~{read2}
    # sort and index
    samtools sort -O bam -o ~{samplename}.sorted.bam ~{samplename}.bam
    samtools index ~{samplename}.sorted.bam
  >>>

  output {
    File sorted_bam = "~{samplename}.sorted.bam"
    File sorted_bam_bai = "~{samplename}.sorted.bam.bai"
  }

  runtime {
    docker: "staphb/smalt:0.7.6"
    memory: "2 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 1
  }
}