# germcall
```
usage: germcall.py [-h] [-b BAMS [BAMS ...]] [-c CONFIG] [-e EXOMEBED]
                   [-l LIST_OF_TOOLS [LIST_OF_TOOLS ...]] [-o OUTPUT_DIR]
                   [-d DRYRUN]

Run germline variation calling using various pipelines

optional arguments:
  -h, --help            show this help message and exit
  -b BAMS [BAMS ...], --bams BAMS [BAMS ...]
                        Bam files to run germline calls
  -c CONFIG, --config CONFIG
  -e EXOMEBED, --exomeBed EXOMEBED
                        Exome capture BED file, supply only if Bam files
                        specified in -b arguments are whole exome sequencing
  -l LIST_OF_TOOLS [LIST_OF_TOOLS ...], --list_of_tools LIST_OF_TOOLS [LIST_OF_TOOLS ...]
                        List of tools to run, by default runs all available
                        tools. Curently supported tools: strelka, freebayes,
                        gatk. If running multiple tools, simply append
                        additional tool name separated with a space. ex) -l
                        strelka gatk
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        output directory to write results
  -d DRYRUN, --dryrun DRYRUN
                        If set to 1, will print out shell commands that ar
                        run. If set to 0 will actually run the commands
```

