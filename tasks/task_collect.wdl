version 1.0

task collect_files {
  input{
    Array[File] bams
    Array[File] bais
    Array[File] vcfs
  }

  command <<<
    mkdir -p intermediate_files
    for f in ~{sep=" " bams}; do
    cp $f intermediate_files/
    done
    for f in ~{sep=" " bais}; do
    cp $f intermediate_files/
    done
    for f in ~{sep=" " vcfs}; do
    cp $f intermediate_files/
    done

    tar -czf intermediate_files.tar.gz intermediate_files
  >>>

  output {
    File intermediate_files = "intermediate_files.tar.gz"
  }

  runtime {
    docker: "ubuntu:jammy-20240911.1"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 1
  }
}