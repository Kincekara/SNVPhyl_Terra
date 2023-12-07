version 1.0

task index {
  input {
    File genome
    String filename = basename(genome)    
  }

  runtime {
    docker: "staphb/samtools:1.9"
    memory: "2 GB"
    cpu: 1
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 1
  } 
}