version 1.0

task freebayes {
  input {
    String samplename
    File sorted_bam
    File reference
    File reference_fai
  }

  command <<<
    freebayes --bam ~{sorted_bam} --ploidy 1 --fasta-reference ~{reference} --vcf ~{samplename}.freebayes.vcf
  >>>

  output {
    File freebayes_vcf = "~{samplename}.freebayes.vcf"
  }

  runtime {
    docker: "staphb/freebayes:1.3.6"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 1
  }
}

