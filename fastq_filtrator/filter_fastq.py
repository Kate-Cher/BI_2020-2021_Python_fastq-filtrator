import sys


def length_filter(minimal, keep):
    # some code here
    print("min_len")


# Parsing arguments part

keep_filtered = 0
max_gc = 100
get_output_name = 0
i = 1
while i < len(sys.argv):
    if sys.argv[i] == "--min_length":
        min_length = sys.argv[i + 1]
        i += 1
    if sys.argv[i] == "--keep_filtered":
        keep_filtered = 1
    if sys.argv[i] == "--gc_bounds":
        min_gc = sys.argv[i + 1]
        i += 1
        if int(sys.argv[i + 1]):
            max_gc = sys.argv[i + 1]
            i += 1
    if sys.argv[i] == "--output_base_name":
        output_base_name = sys.argv[i + 1]
        i += 1
        get_output_name = 1
    if sys.argv[i][-6:] == '.fastq':
        input_file = sys.argv[i]
    i += 1

output_base_name_passed = output_base_name + "__passed.fastq"
output_base_name_failed = output_base_name + "--failed.fastq"




