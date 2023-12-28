version 1.0

import "../tasks/task_smalt.wdl" as smalt
import "../tasks/task_snvphyl_tools.wdl" as tools
import "../tasks/task_vcf2snv.wdl" as vcf2snv
import "../tasks/task_phyml.wdl" as phyml
import "../tasks/task_report.wdl" as report

import "wf_variants.wdl" as variants

workflow snvphyl {
  input {
    File inputSamplesFile
    Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
    File reference
    Int window_size = 11
    Int density_threshold = 2
    String? colorscale
    Int? tree_width
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

  call tools.verify_map_q {
    input:
      sorted_bams = variants.sorted_bam,
      sorted_bam_bais = variants.sorted_bam_bai
  }

  call vcf2snv.vcf2snv {
    input:
      filtered_densities = variants.filtered_density,
      consolidated_bcfs = variants.consolidated_bcf,
      consolidated_bcfs_csis = variants.consolidated_bcf_csi,
      invalid_positions = find_repeats.invalid_positions,
      reference = index.reference
  }

  call tools.stats_and_matrix {
    input:
      snvtable = vcf2snv.snvtable,
      snvalignment = vcf2snv.snvalignment
  }

  call phyml.phyml {
    input:
      snvalignment = vcf2snv.snvalignment
  }

  call report.create_report {
    input:
      matrix = stats_and_matrix.snvmatrix,
      reference = reference,
      vcf2core = vcf2snv.vcf2core,
      newick = phyml.tree,
      colorscale = colorscale,
      tree_width = tree_width
  }

  output {
    File mapping_quality = verify_map_q.mapping_quality
    File vcf2core = vcf2snv.vcf2core
    File filter_stats = stats_and_matrix.filterstats
    File snv_matrix = stats_and_matrix.snvmatrix
    File phyml_tree = phyml.tree
    File phyml_tree_stats = phyml.treestats
    File summary_report = create_report.summary_report
  }
}
