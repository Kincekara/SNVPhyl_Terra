version 1.0

task vcf2snv {
  input {
    Array[File] filtered_densities
    Array[File] consolidated_bcfs
    Array[File] consolidated_bcfs_csis
    File invalid_positions
    File reference
  }

  command <<<
    # merge filtered densities
    fpaths=(~{sep=' ' filtered_densities})
    for file in ${fpaths[@]}; do
      cat $file
    done > filtered_density_all.txt

    # merge filtered densities and invalid positions
    cat filtered_density_all.txt ~{invalid_positions}> new_invalid_positions.bed

    # prep bcfs
    bpath=(~{sep=' ' consolidated_bcfs})
    for f in ${bpath[@]}; do
      fname=$(basename $f .consolidated.bcf)
      echo "--consolidate_vcf $fname=$f " | tr -d "\n" >> consolidation_line.txt
    done
    consolidate_cmd=$(cat consolidation_line.txt)

    # vcf2snv
    vcf2snv_alignment.pl --reference reference --invalid-pos new_invalid_positions.bed --format fasta --format phylip --numcpus 4 --output-base snvalign --fasta ~{reference} $consolidate_cmd
    mv snvalign-positions.tsv snvTable.tsv
    mv snvalign-stats.csv vcf2core.tsv
    if [[ -f snvalign.phy ]]; then
        mv snvalign.phy snvAlignment.phy
        sed -i "s/'//" snvAlignment.phy
        sed -i "s/'//" snvAlignment.phy
    else
        touch snvAlignment.phy
    fi  
  >>>

  output {
    File snvtable = "snvTable.tsv"
    File vcf2core = "vcf2core.tsv"
    File snvalignment = "snvAlignment.phy"
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

