version 1.0

import "../tasks/task_smalt.wdl" as smalt
import "../tasks/task_snvphly_tools.wdl" as tools
import "../tasks/task_freebayes.wdl" as freebayes
import "../tasks/task_bcftools.wdl" as bcftools

workflow variants {

  input {
    File reference
    String samplename
    File read1
    File read2
    Int window_size
    Int density_threshold
    File genome
    File fai
    File sma
    File smi
  }

  call smalt.map {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2,
      genome = genome,
      fai = fai,
      sma = sma,
      smi = smi
  }
  
  call freebayes.freebayes {
    input:
      samplename = samplename,
      sorted_bam = map.sorted_bam,
      reference = reference
  }

  call tools.filter_and_bcf {
    input:
      samplename = samplename,
      freebayes_vcf = freebayes.freebayes_vcf
  }

  call bcftools.mpileup {
    input:
      samplename = samplename,
      sorted_bam = map.sorted_bam,
      reference = reference
  }

  call tools.consolidate_bcf {
    input:
      samplename = samplename,
      mpileup_bcf = mpileup.mpileup_bcf,
      filtered_bcf = filter_and_bcf.filtered_bcf,
      filtered_csi = filter_and_bcf.filtered_csi,
      window_size = window_size,
      density_threshold = density_threshold
  }
 
  output {
    File sorted_bam = map.sorted_bam
    File sorted_bam_bai = map.sorted_bam_bai
    File consolidated_bcf = consolidate_bcf.consolidated_bcf
    File consolidated_vcf = consolidate_bcf.consolidated_vcf
    File consolidated_bcf_csi = consolidate_bcf.consolidated_bcf_csi
    File filtered_density = consolidate_bcf.filtered_density
  }
}
