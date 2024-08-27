version 1.0

task create_report {
  input {
    File reference
    File vcf2core
    File matrix
    File newick
    String colorscale = "YlGnBu_r"
    Int tree_width = 800
  }

  command <<<
    python3 <<CODE
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import base64
    import toytree      
    import toyplot       
    import re     
    
    # read matrix
    df = pd.read_csv("~{matrix}", sep = '\t', index_col="strain")
    df.dropna(how='all', axis=1, inplace=True)
    df2 = df.drop(columns="reference", index="reference")
    df2 = df2.rename_axis(None, axis=0)

    # read vcf2core
    core_df = pd.read_csv("~{vcf2core}", sep = '\t', index_col="#Reference name")
    core_percent = core_df.loc["all"]["Percentage of valid and included positions in core genome"]

    # read fasta header
    with open("~{reference}", "r") as f:
        fasta_header = f.readlines()[0].replace(">", "").strip()

    # plot matrix
    dim = df2.shape[0] * 0.6
    fig, ax = plt.subplots(figsize=(dim, dim))
    ax = sns.heatmap(df2, cmap="~{colorscale}", annot=True, cbar=False, fmt="g", linewidths=0.5, square=True)
    ax.xaxis.tick_top()
    plt.xticks(rotation=30, ha="left", rotation_mode="anchor")
    plt.savefig('heatmap.png', bbox_inches='tight')
    with open('heatmap.png','rb') as f:
        b64data = base64.b64encode(f.read())
    bstring = b64data.decode()

    ## draw phylogenetic tree
    with open("~{newick}", "r") as f:
        nwk = f.read()

    # remove reference from tree
    nwk = re.sub('reference:[0-9]+\\.[0-9]+,', '', nwk)
    nwk = re.sub(',reference:[0-9]+\\.[0-9]+', '', nwk)

    # draw tree
    tree = toytree.tree(nwk)
    canvas, axes, mark = tree.draw(tip_labels_style={"font-size": "13px"}, scale_bar=True, width=~{tree_width})
    toyplot.html.render(canvas, "tree-plot.html")

    # report style
    style=''' <style>
            h1 {
              font-family: Arial, Helvetica, sans-serif;
              font-size: 1.5em;
              text-align: center;
              color: #1e5c85
            }

            h2 {
              font-family: Arial, Helvetica, sans-serif;
              font-size: 1em;
              margin-top: 1em;
              margin-bottom: 0.4em;
              text-align: left;
              color: #D6672E 
            }

            table {
              font-family: Arial, Helvetica, sans-serif;
              font-size: 0.8em;
              border: 1px solid;
            }
            
            th {
              background-color: #f1f1f1;
            }

            th, td {
              text-align: left;
              padding-left: 5px;
              padding-right: 5px;
              padding-top: 5px;
              padding-bottom: 5px;
            }        

            img {
              max-width: 100%;
              height: auto;
            }

            footer {
              font-family: Arial, Helvetica, sans-serif;
              font-size: 0.8em;
            }   
            </style>'''

    footer = '<p><i>This report is created by <a href="https://github.com/Kincekara/SNVPhyl_Terra">SNVPhyl_Terra</a> bioinformatics pipeline.</br></i></p>'

    with open("tree-plot.html", "r") as f:
        html_tree = f.read()

    # generate report
    with open("snvphyl_report.html", "w") as f:
        f.write('<!DOCTYPE html><html lang="en-us"><head><meta charset="UTF-8">\n<title>SNVPhyl report</title>\n')
        f.write(style + "\n</head>\n")
        f.write("<body>\n<header><h1>SNVPhyl Report</h1></header>\n<hr>\n<article>\n")
        f.write("<h2>Stats</h2>\n")
        f.write('<p><table><th></th><th></th><tr><td>Reference</td><td>' + fasta_header + "</td><tr>")
        f.write("<tr><td>SNVPhyl core genome</td><td>" + str(core_percent) + "%</td></tr></table><p>\n")
        f.write("<h2>SNV Matrix</h2>\n")
        f.write('<img src="data:image/png;base64,' + bstring + '"' + 'alt="Heatmap" />\n')
        f.write("<h2>Phylogenetic Tree</h2>\n")
        f.write("<p><p>\n")
        f.write(html_tree)
        f.write("\n</article>\n<hr>\n<footer>" + footer + "</footer>\n</body>\n</html>\n")
    CODE
  >>>

  output {
    File summary_report = "snvphyl_report.html"
  }

  runtime {
    docker: "kincekara/python-tools:0.2"
    memory: "2 GB"
    cpu: 1
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 1
  }
}