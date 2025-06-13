version 1.0

task map {
  input {
    String samplename
    File read1
    File read2    
    File genome
    String filename = basename(genome)        
    Int min_insert_size = 20
    Int max_insert_size = 1000
  }

  command <<<  
    # copy reference to local directory
    cp ~{genome} ./~{filename}
    # index reference
    smalt index -k 13 -s 6 ~{filename} ~{genome}    
    samtools faidx ~{filename}  
    # map reads
    smalt map -f bam -n 4 -l pe -i ~{max_insert_size} -j ~{min_insert_size} -r 1 -y 0.5 -o ~{samplename}.bam ~{filename} ~{read1} ~{read2}
    # sort and index
    samtools sort -O bam -o ~{samplename}.sorted.bam ~{samplename}.bam
    samtools index ~{samplename}.sorted.bam
  >>>

  output {
    File sorted_bam = "~{samplename}.sorted.bam"
    File sorted_bam_bai = "~{samplename}.sorted.bam.bai"
    File reference = "~{filename}"
    File reference_fai = "~{filename}.fai"
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