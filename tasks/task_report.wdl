version 1.0

task create_report {
  input {
    File reference
    File vcf2core
    File matrix
    File newick
    String colorscale = "YlGn_r"
    Int tree_width = 600
  }

  command <<<
    python3 <<CODE
    import pandas as pd
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
        fasta_header = f.readlines()[0].replace(">", "")

    # plot matrix
    styled_table = df2.style.background_gradient(cmap="~{colorscale}", axis=None)
    styled_table.set_table_styles([{'selector': 'td', 'props': 'text-align: center;'}])
    html_table = styled_table.to_html()

    ## draw phylogenetic tree
    with open("~{newick}", "r") as f:
        nwk = f.read()

    # remove reference from tree
    nwk = re.sub('reference:[0-9]+\\.[0-9]+,', '', nwk)
    nwk = re.sub(',reference:[0-9]+\\.[0-9]+', '', nwk)

    # draw tree
    tre = toytree.tree(nwk, tree_format=0)
    canvas, axes, mark = tre.draw(tip_labels_style={"font-size": "13px"}, scalebar=True, width=~{tree_width})
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
        f.write('<!DOCTYPE html><html lang="en-us"><head><meta charset="UTF-8"><title>SNVPhyl report</title></head>')
        f.write(style + "</head>")
        f.write("<body><header><h1>SNVPhyl Report</h1></header><hr><article>")
        f.write("<h2>Stats</h2>")
        f.write('<p><table><th></th><th></th><tr><td>Reference</td><td>' + fasta_header + "</td><tr>")
        f.write("<tr><td>SNVPhyl core genome</td><td>" + str(core_percent) + "%</td></tr></table><p>")
        f.write("<h2>SNV Matrix</h2>")
        f.write(html_table)
        f.write("<h2>Phylogenetic Tree</h2>")
        f.write("<p><p>")
        f.write(html_tree)
        f.write("</article><hr><footer>" + footer + "</footer></body></html>")
    CODE
  >>>

  output {
    File summary_report = "snvphyl_report.html"
  }

  runtime {
    docker: "kincekara/python-tools:0.1"
    memory: "2 GB"
    cpu: 1
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 1
  }
}