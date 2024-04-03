version 1.0

import "../tasks/task_prep.wdl" as prep
import "wf_snvphyl_local.wdl" as main

workflow snvphyl_terra {

  meta {
    description: "Terra implementation of SNVPhyl pipeline"  
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

  if ( validate_inputs.check == "PASS" ) {
    call main.snvphyl {
      input:
      inputSamplesFile = validate_inputs.samplelist,
      reference = validate_inputs.ref_genome,
      window_size = window_size,
      density_threshold = density_threshold,
      colorscale = colorscale,
      tree_width = tree_width
    }
  }

  output {
    String version = "SNVPhyl_Terra v1.0.3-dev"
    String inputs_check = validate_inputs.check
    File? mapping_quality = snvphyl.mapping_quality
    File? vcf2core = snvphyl.vcf2core
    File? filter_stats = snvphyl.filter_stats
    File? snv_matrix = snvphyl.snv_matrix
    File? phyml_tree = snvphyl.phyml_tree
    File? phyml_tree_stats = snvphyl.phyml_tree_stats
    File? summary_report = snvphyl.summary_report
  }
}