version 1.0

import "../tasks/task_smalt.wdl" as smalt
import "../tasks/task_snvphly_tools.wdl" as tools
import "../tasks/task_freebayes.wdl" as freebayes
import "../tasks/task_bcftools.wdl" as bcftools

workflow snvphyl {

  meta {
	description: "SNVPhyl"
  }

  input {
    File reference
    String samplename
    File read1
    File read2
    Int window_size = 11
    Int density_threshold = 2
  }

  call smalt.index {
    input:
      genome = reference
  }

  call tools.find_repeats {
    input:
      genome = reference
  }

  call smalt.map {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2,
      genome = index.reference,
      fai = index.fai,
      sma = index.sma,
      smi = index.smi
  }
  
  call tools.verify_map_q {
    input:
      samplename = samplename,
      sorted_bam = map.sorted_bam,
      sorted_bam_bai = map.sorted_bam_bai
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
  }


}
