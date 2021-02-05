import sys
import os.path


def length_filter(read, minimal):
    if len(read) >= int(minimal):
        return 1
    else:
        return 0


def gc_bounder(read, gc_min, gc_max):
    gc_count = (read.count('C') + read.count('G') + read.count('c') + read.count('g')) / (len(read) - 1) * 100
    if gc_min <= gc_count <= gc_max:
        return 1
    else:
        return 0


# Optional part with help message (maybe useful)

if sys.argv[1] == '-h' or sys.argv[1] == '--help':
    print("This is fastq reads filtrator. \n"
          "Parameters you can use: \n"
          "\t --min_length:         Minimal length of read to keep (Optional argument)\n"
          "\t --keep_filtered:      Use this optional argument to save filtered reads \n"
          "\t --output_base_name:   Use this optional argument to set the base of the output filename \n"
          "\t\t\t\t'--output_base_name out_file' output files called out_file__passed.fastq (+ out_file__failed) \n"
          "\t --gc_bounds:          Use this optional arg to set GC content \n"
          "\t\t\t\t'--gc_bounds 55'     saves reads with GC content >= 55%\n"
          "\t\t\t\t'--gc_content 55 70' saves reads with GC content >=55% & <=70%\n"
          "fastq file is required positional argument\n\n"
          "Using example:\n"
          "\tpython filter_fastq.py --some_argument some_values test.fastq")
else:
    keep_filtered = 0
    max_gc = 100.0
    min_gc = 0.0
    min_length = 0
    i = 1
    type_er_fl = 0
    prohib_symb = ['?', '/', ':', '*', '>', '<', '|', '+', '"', '\\', '--min_length', '--keep_filtered',
                   '--gc_bounds']

    if sys.argv[len(sys.argv) - 1][-6:] == '.fastq':
        input_file = sys.argv[len(sys.argv) - 1]
        output_base_name = input_file[:-6]

# Parsing arguments
        if os.path.isfile(input_file):

            while i < len(sys.argv):
                if sys.argv[i] == "--min_length":
                    min_length = sys.argv[i + 1]
                    i += 1
                if sys.argv[i] == "--keep_filtered":
                    keep_filtered = 1
                if sys.argv[i] == "--gc_bounds":
                    try:
                        min_gc = float(sys.argv[i + 1])
                    except ValueError:
                        type_er_fl = 1
                    try:
                        max_gc = float(sys.argv[i + 2])
                        i += 1
                    except ValueError:
                        pass
                    i += 1
                if sys.argv[i] == "--output_base_name":
                    f = 1
                    for symb in prohib_symb:
                        if symb in sys.argv[i + 1]:
                            f = 0
                    if f == 1:
                        output_base_name = sys.argv[i + 1]
                        i += 1
                    else:
                        print("Wrong base name. Using default value \n"
                              f"Prohibited symbols: ")
                        print(*prohib_symb)
                        print('Using default base')
                i += 1

            output_base_name_passed = output_base_name + "__passed.fastq"
            output_base_name_failed = output_base_name + "__failed.fastq"

            # Run information module. Checking argument correctness

            print('Run information: ')
            print(f"\tPassed reads are in {output_base_name_passed}")
            if keep_filtered == 1:
                print(f"\tFailed reads are in {output_base_name_failed}")
            else:
                print('\tReads that have not passed the filter are not saved')
            print(f"Passed reads: \n"
                  f"\t Minimal read's length: {min_length}")
            if type_er_fl == 1:
                print("Error in GC content type. Must be float. Using default values")
            else:
                if max_gc < min_gc <= 100:
                    print('Minimal GC content must be first. Changing positions')
                    min_gc, max_gc = max_gc, min_gc
                if max_gc > 100:
                    print('There were error in maximal GC content (> 100%). Using default value')
                    max_gc = 100
                if min_gc > 100:
                    print('There were error in minimal GC content (> 100%). Using default value')
                    min_gc = 0
            print(f"\t Minimal read's gc content: {min_gc} \n"
                  f"\t Maximal read's gc content: {max_gc}\n")

            # Checking reads part

            passed_counter = 0
            failed_counter = 0
            with open(input_file) as fastq:
                with open(output_base_name_passed, 'w') as output_passed:
                    with open(output_base_name_failed, 'w') as output_failed:
                        for line in fastq:
                            lines = [line]
                            for i in range(3):
                                lines.append(fastq.readline())
                            if length_filter(lines[1], min_length) * gc_bounder(lines[1], min_gc, max_gc) == 1:
                                output_passed.writelines(lines)
                                passed_counter += 1
                            else:
                                failed_counter += 1
                                if keep_filtered == 1:
                                    output_failed.writelines(lines)

            all_reads = passed_counter + failed_counter
            print(f"Number of passed filtration reads: {passed_counter}")
            print("Percent of passed filtration reads: %.4f" % (passed_counter / all_reads * 100), "%")
            print(f"\nNumber of failed filtration reads: {failed_counter}")
            print("Percent of failed filtration reads: %.4f" % (failed_counter / all_reads * 100), "%")
        else:
            print('No such file in directory')
    else:
        print("There is no fastq file. Check file name correctness and try again. \n"
              "Fastq file is required last positional argument")
