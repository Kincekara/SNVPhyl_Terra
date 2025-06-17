version 1.0

import "../tasks/task_snvphyl_tools.wdl" as tools
import "../tasks/task_vcf2snv.wdl" as vcf_to_snv
import "../tasks/task_phyml.wdl" as tree
import "../tasks/task_report.wdl" as report
import "../tasks/task_collect.wdl" as collect
import "wf_variants.wdl" as snv

workflow snvphyl {
  input {
    File inputSamplesFile
    Array[Array[String]] inputSamples = read_tsv(inputSamplesFile)
    File reference
    Int window_size = 11
    Int density_threshold = 2
    Int min_coverage = 10
    Int min_mean_mapping = 30
    Float snv_abundance_ratio = 0.75
    String? colorscale
    Int? tree_width
    Boolean intermediate_files = false
  }

  call tools.find_repeats {
    input:
      genome = reference
  }

  scatter (sample in inputSamples) {
    call snv.variants {
      input:
        samplename = sample[0],
        read1 = sample[1],
        read2 = sample[2],
        reference = reference,
        window_size = window_size,
        density_threshold = density_threshold,
        min_coverage = min_coverage,
        min_mean_mapping = min_mean_mapping,
        snv_abundance_ratio = snv_abundance_ratio
    }
  }

  call tools.verify_map_q {
    input:
      sorted_bams = variants.sorted_bam,
      sorted_bam_bais = variants.sorted_bam_bai,
      min_depth = min_coverage
  }

  call vcf_to_snv.vcf2snv {
    input:
      filtered_densities = variants.filtered_density,
      consolidated_bcfs = variants.consolidated_bcf,
      consolidated_bcfs_csis = variants.consolidated_bcf_csi,
      invalid_positions = find_repeats.invalid_positions,
      reference = reference
  }

  call tools.stats_and_matrix {
    input:
      snvtable = vcf2snv.snvtable,
      snvalignment = vcf2snv.snvalignment
  }

  call tree.phyml {
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

  if (intermediate_files) {
    call collect.collect_files {
      input:
        bams = variants.sorted_bam,
        bais = variants.sorted_bam_bai,
        vcfs = variants.consolidated_vcf
    }
  }

  output {
    String version = "SNVPhyl_Terra v1.0.7"
    File mapping_quality = verify_map_q.mapping_quality
    File vcf2core = vcf2snv.vcf2core
    File filter_stats = stats_and_matrix.filterstats
    File snv_alignment = vcf2snv.snvalignment
    File snv_matrix = stats_and_matrix.snvmatrix
    File phyml_tree = phyml.tree
    File phyml_tree_stats = phyml.treestats
    File summary_report = create_report.summary_report
    File? bams_and_vcfs = collect_files.intermediate_files
  }
}
