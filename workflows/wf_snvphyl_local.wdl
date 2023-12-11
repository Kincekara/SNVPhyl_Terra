version 1.0

import "../tasks/task_smalt.wdl" as smalt
import "../tasks/task_snvphly_tools.wdl" as tools
import "../tasks/task_vcf2snv.wdl" as vcf2snv

import "wf_variants.wdl" as variants

workflow snvphyl {
  input {
    File inputSamplesFile
    Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
    File reference
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

  scatter (sample in inputSamples) {
    call variants.variants {
      input:
        samplename = sample[0],
        read1 = sample[1],
        read2 = sample[2],
        reference = reference,
        genome = index.reference,
        fai = index.fai,
        sma = index.sma,
        smi = index.smi,        
        window_size = window_size,
        density_threshold = density_threshold

    }
  }

  call vcf2snv.vcf2snv {
    input:
      filtered_densities = variants.filtered_density,
      consolidated_bcfs = variants.consolidated_bcf,
      invalid_positions = find_repeats.invalid_positions,
      reference = index.reference
  }
}
