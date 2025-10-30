version 1.0

task phastest {
  input {
    File genome
  }

command <<<
  # find prophage sequences
  phastest -i contig -s ~{genome} --silent --phage-only --yes

  # prep outputs
  echo -e "#Chromosome\tStart\tEnd" > phage_sites.bed
  jq -r '.[] | select(.completeness == "intact") | [(.contig_tag | split(",")[0]), .start, .stop] | @tsv' ./phastest_out/*/json_input_regions >> phage_sites.bed
  cp ./phastest_out/*/json_input_regions ./phage_regions.json
>>>

  output {
    File phage_bed = "phage_sites.bed"
    File phage_json = "phage_regions.json"
  }

  runtime {
    docker: "kincekara/phastest-cli:0.1"
    memory: "16 GB"
    cpu: 8
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 1
  }
}