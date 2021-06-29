# BI_2020-2021_Python_fastq-filtrator
Fastq filtering "tool".

fastq-filtrator is a small console tool for filtering fastq reads.

Parameters you can use:
- --min_length:         Minimal length of read to keep (Optional argument)
- --keep_filtered:      Use this optional argument to save filtered reads 
- --output_base_name:   Use this optional argument to set the base of the output filename
- --output_base_name out_file' output files called out_file__passed.fastq (+ out_file__failed)
    - --gc_bounds:          Use this optional arg to set GC content
    - --gc_bounds 55'     saves reads with GC content >= 55%
    - --gc_content 55 70' saves reads with GC content >=55% & <=70%
- fastq file is required positional argument
