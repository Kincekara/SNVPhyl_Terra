version 1.0

task phyml {
  input {
    File snvalignment
    String filename = basename(snvalignment) 
  }

  command <<<
    cp ~{snvalignment} ./~{filename}
    phyml -i ~{filename} --datatype nt --model GTR -v 0.0 -s BEST --ts/tv e --nclasses 4 --alpha e --bootstrap -4 --quiet
    mv ~{filename}_phyml_stats.txt phylogeneticTreeStats.txt
    mv ~{filename}_phyml_tree.txt phylogeneticTree.newick
  >>>

  output {
    File treestats = "phylogeneticTreeStats.txt"
    File tree = "phylogeneticTree.newick"
  }

  runtime {
    docker: "staphb/phyml:3.3.20220408"
    memory: "2 GB"
    cpu: 1
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 1
  }
}