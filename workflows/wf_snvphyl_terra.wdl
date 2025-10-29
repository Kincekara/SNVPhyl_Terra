version 1.0

import "../tasks/task_prep.wdl" as prep
import "wf_snvphyl_local.wdl" as main

workflow snvphyl_terra {

  meta {
    author: "Kutluhan Incekara"
    email: "kutluhan.incekara@ct.gov"
    description: "Terra implementation of SNVPhyl pipeline"  
  }

  parameter_meta {
    samplename: {
      description: "Array of sample names"
    }
    read1: {
      description: "Array of read1 files"
    }
    read2: {
      description: "Array of read2 files"
    }
    taxon: {
      description: "Taxon name"
    }
    reference: {
      description: "Reference genome file",
      optional: "true",
      patterns: ["*.fasta", "*.fa", "*.fna"]
    }
    accession: {
      description: "NCBI Accession number for reference genome",
      optional: "true"
    }
    window_size: {
      description: "Window size for SNVPhyl",
      optional: "true",
      default: "11"
    }
    density_threshold: {
      description: "Density threshold for SNVPhyl",
      optional: "true",
      default: "2"
    }
    colorscale: {
      description: "Color scale for SNVPhyl",
      optional: "true",
      default: "YlGnBu_r"
    }
    tree_width: {
      description: "Tree width for SNVPhyl",
      optional: "true",
      default: "800"
    }
    intermediate_files: {
      description: "Whether to collect bam and vcf files into a tarball",
      optional: "true",
      default: "false"
    }
  }

  input {
    Array[String] samplename
    Array[File] read1
    Array[File] read2    
    String taxon
    File? reference
    String? accession
    Int? window_size
    Int? density_threshold
    String? colorscale
    Int? tree_width
    Boolean? intermediate_files
  }

  call prep.validate_inputs {
    input:
    samplename = samplename,
    read1 = read1,
    read2 =read2,
    reference = reference,
    taxon = taxon,
    accession = accession
  }

  call main.snvphyl {
    input:
    inputSamplesFile = validate_inputs.samplelist,
    reference = validate_inputs.ref_genome,
    window_size = window_size,
    density_threshold = density_threshold,
    colorscale = colorscale,
    tree_width = tree_width,
    intermediate_files = intermediate_files,
  }

  output {
    String version = snvphyl.version
    String inputs_check = validate_inputs.check
    File mapping_quality = snvphyl.mapping_quality
    File vcf2core = snvphyl.vcf2core
    File filter_stats = snvphyl.filter_stats
    File snv_matrix = snvphyl.snv_matrix
    File snv_alignment = snvphyl.snv_alignment
    File phyml_tree = snvphyl.phyml_tree
    File phyml_tree_stats = snvphyl.phyml_tree_stats
    File summary_report = snvphyl.summary_report
    File? bams_and_vcfs = snvphyl.bams_and_vcfs
  }
}