version 1.0

task vcf_to_bcf {
  input{
    String samplename
    File filtered_vcf
  }

  command <<<
    bcftools index -f ~{filtered_vcf}
    bcftools view --output-type b --output-file ~{samplename}.freebayes.filtered.bcf ~{filtered_vcf}
    bcftools index -f ~{samplename}freebayes_filtered.bcf 
  >>>

  output {
    File filtered_bcf = "~{samplename}freebayes_filtered.bcf"
    File filtered_csi = "~{samplename}freebayes_filtered.bcf.csi"
  }

  runtime {
    docker: "staphb/bcftools:1.15"
    memory: "1 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 1
  }
}

task mpileup {
  input {
    String samplename
    File sorted_bam
    File reference
  }

  command <<<
    bcftools mpileup --threads 4 --fasta-ref ~{reference} -A -B -C 0 -d 1024 -q 0 -Q 0 --output-type b -I --output ~{samplename}.mpileup.vcf.gz ~{sorted_bam}
    bcftools index -f ~{samplename}.mpileup.vcf.gz
    bcftools call --ploidy 1 --threads 4 --output ~{samplename}.mpileup.bcf --output-type b --consensus-caller ~{samplename}.mpileup.vcf.gz
  >>>

  output {
    File mpileup_bcf = "~{samplename}.mpileup.bcf"
  }

  runtime {
    docker: "staphb/bcftools:1.15"
    memory: "1 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 1
  }
}