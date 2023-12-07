version 1.0

task find_repeats {
  input {
    File genome
  }

  command <<<
    find-repeats.pl ~{genome} --min-length 150 --min-pid 90 > invalid_positions.bed
  >>>

  output {
    File invalid_positions = "invalid_positions.bed"
  }

  runtime {
    docker: "staphb/snvphyl-tools:1.8.2"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 1
  }
}

task verify_map_q {
  input {
    String samplename
    File sorted_bam
    File sorted_bam_bai
  }

  command <<<    
    verify_mapping_quality.pl -c 4 --min-depth 10 --min-map 80 --output ~{samplename}.mappingQuality.txt --bam bam1=~{sorted_bam}
  >>>

  output {
    File mapping_quality = "~{samplename}.mappingQuality.txt"
  }

  runtime {
    docker: "staphb/snvphyl-tools:1.8.2"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 1
  }

}

task filter_and_bcf {
  input {
    String samplename
    File freebayes_vcf
  }

  command <<<
    # filter
    filterVcf.pl --noindels ~{freebayes_vcf} -o ~{samplename}.freebayes.filtered.vcf
    # compress 
    bgzip -f ~{samplename}.freebayes.filtered.vcf
    # vcf to bcf
    bcftools index -f ~{samplename}.freebayes.filtered.vcf.gz
    bcftools view --output-type b --output-file ~{samplename}.freebayes.filtered.bcf ~{samplename}.freebayes.filtered.vcf.gz
    bcftools index -f ~{samplename}.freebayes.filtered.bcf 

  >>>

  output {
    File filtered_vcf = "~{samplename}.freebayes.filtered.vcf.gz"
    File filtered_bcf = "~{samplename}.freebayes.filtered.bcf"
    File filtered_csi = "~{samplename}.freebayes.filtered.bcf.csi"
  }

  runtime {
    docker: "staphb/snvphyl-tools:1.8.2"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 1
  }
}

task consolidate_bcf {
  input {
    String samplename
    File mpileup_bcf
    File filtered_bcf
    File filtered_csi
    Int window_size
    Int density_threshold
  }

  command <<<
    consolidate_vcfs.pl --coverage-cutoff 10 --min-mean-mapping 30 --snv-abundance-ratio 0.75 --vcfsplit ~{filtered_bcf} --mpileup ~{mpileup_bcf} --filtered-density-out ~{samplename}.filtered-density.txt --window-size ~{window_size} --density-threshold ~{density_threshold} -o ~{samplename}.consolidated.bcf > ~{samplename}.consolidated.vcf
    bcftools index -f ~{samplename}.consolidated.bcf
  >>>

  output {
    File consolidated_bcf = "~{samplename}.consolidated.bcf"
    File consolidated_vcf = "~{samplename}.consolidated.vcf"
    File consolidated_bcf_csi = "~{samplename}.consolidated.bcf.csi"
    File filtered_density = "~{samplename}.filtered-density.txt"
  }

  runtime {
    docker: "staphb/snvphyl-tools:1.8.2"
    memory: "30 GB"
    cpu: 1
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 1
  }
}