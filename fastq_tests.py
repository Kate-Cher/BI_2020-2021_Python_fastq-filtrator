import unittest
import os
import sys
from io import StringIO
from copy import deepcopy
from fastq_updated_filter import arg_pars, filter_reads, gc_bounder, length_filter


class TestArgumentParser(unittest.TestCase):
    def setUp(self):
        self.held, sys.stdout = sys.stdout, StringIO()

    def test_default_args(self):
        test_dict = {'keep_filtered': 0, 'max_gc': 100.0, 'min_gc': 0.0,
                     'min_length': 0, 'type_er_fl': 0, 'error_state': False,
                     'output_base_name': 'test', 'input_file': 'test.fastq',
                     'output_base_name_passed': 'test__passed.fastq',
                     'output_base_name_failed': 'test__failed.fastq'}
        arg = ["fastq_updated_filter.py", "test.fastq"]
        self.assertDictEqual(test_dict, arg_pars(arg))

    # Input file name tests
    def test_wrong_file_name(self):
        test_dict = {'keep_filtered': 0, 'max_gc': 100.0, 'min_gc': 0.0,
                     'min_length': 0, 'type_er_fl': 0, 'error_state': True,
                     'output_base_name': '', 'input_file': '',
                     'output_base_name_passed': '',
                     'output_base_name_failed': ''}
        arg = ["fastq_updated_filter.py", "test"]
        self.assertDictEqual(test_dict, arg_pars(arg))
        self.assertEqual("There is no fastq file. Check file name correctness and try again. \n"
                         "Fastq file is required last positional argument\n", sys.stdout.getvalue())

    def test_no_file_arg(self):
        arg = ["fastq_updated_filter.py"]
        arg_pars(arg)
        self.assertEqual("No arguments found. Try -h/--help to help\n", sys.stdout.getvalue())

    def test_no_such_fastq_file(self):
        test_dict = {'keep_filtered': 0, 'max_gc': 100.0, 'min_gc': 0.0,
                     'min_length': 0, 'type_er_fl': 0, 'error_state': True,
                     'output_base_name': '', 'input_file': '',
                     'output_base_name_passed': '',
                     'output_base_name_failed': ''}
        arg = ["fastq_updated_filter.py", "mistake.fastq"]
        arg_pars(arg)
        self.assertEqual("No such file in directory\n", sys.stdout.getvalue())

    # Output file name tests
    def test_wrong_base_names(self):
        # Wrong base name => using default base name (input - ".fastq")
        test_dict = {'keep_filtered': 0, 'max_gc': 100.0, 'min_gc': 0.0,
                     'min_length': 0, 'type_er_fl': 0, 'error_state': False,
                     'output_base_name': 'test', 'input_file': 'test.fastq',
                     'output_base_name_passed': 'test__passed.fastq',
                     'output_base_name_failed': 'test__failed.fastq'}
        arg = ["fastq_updated_filter.py", "--output_base_name", "oof??", "test.fastq"]
        self.assertDictEqual(test_dict, arg_pars(arg))

    def test_no_base(self):
        test_dict = {'keep_filtered': 0, 'max_gc': 100.0, 'min_gc': 0.0,
                     'min_length': 0, 'type_er_fl': 0, 'error_state': False,
                     'output_base_name': 'test', 'input_file': 'test.fastq',
                     'output_base_name_passed': 'test__passed.fastq',
                     'output_base_name_failed': 'test__failed.fastq'}
        arg = ["fastq_updated_filter.py", "test.fastq"]
        self.assertDictEqual(test_dict, arg_pars(arg))

    # Testing parsing correctness
    def test_min_length_arg(self):
        test_dict = {'keep_filtered': 0, 'max_gc': 100.0, 'min_gc': 0.0,
                     'min_length': 50, 'type_er_fl': 0, 'error_state': False,
                     'output_base_name': 'test', 'input_file': 'test.fastq',
                     'output_base_name_passed': 'test__passed.fastq',
                     'output_base_name_failed': 'test__failed.fastq'}
        # Right minimal length
        arg = ["fastq_updated_filter.py", "--min_length", "50", "test.fastq"]
        self.assertDictEqual(test_dict, arg_pars(arg))

        # Mistake in min_length (using default)
        arg = ["fastq_updated_filter.py", "--min_length", "50ooh", "test.fastq"]
        test_dict["min_length"] = 0
        self.assertDictEqual(test_dict, arg_pars(arg))

    def test_keep_filtered(self):
        test_dict = {'keep_filtered': 1, 'max_gc': 100.0, 'min_gc': 0.0,
                     'min_length': 0, 'type_er_fl': 0, 'error_state': False,
                     'output_base_name': 'test', 'input_file': 'test.fastq',
                     'output_base_name_passed': 'test__passed.fastq',
                     'output_base_name_failed': 'test__failed.fastq'}
        arg = ["fastq_updated_filter.py", "--keep_filtered", "test.fastq"]
        self.assertDictEqual(test_dict, arg_pars(arg))

    def test_gc_min(self):
        test_dict = {'keep_filtered': 0, 'max_gc': 100.0, 'min_gc': 40.0,
                     'min_length': 0, 'type_er_fl': 0, 'error_state': False,
                     'output_base_name': 'test', 'input_file': 'test.fastq',
                     'output_base_name_passed': 'test__passed.fastq',
                     'output_base_name_failed': 'test__failed.fastq'}
        # Correct gc min
        arg = ["fastq_updated_filter.py", "--gc_bounds", "40", "test.fastq"]
        self.assertDictEqual(test_dict, arg_pars(arg))

        # Mistake min gc (>100). Using default (0)
        arg = ["fastq_updated_filter.py", "--gc_bounds", "140", "test.fastq"]
        test_dict["min_gc"] = 0
        self.assertDictEqual(test_dict, arg_pars(arg))

        # Mistake in gc type
        arg = ["fastq_updated_filter.py", "--gc_bounds", "--holy(", "test.fastq"]
        test_dict["type_er_fl"] = 1
        self.assertDictEqual(test_dict, arg_pars(arg))

    def test_gc_both_bounds(self):
        test_dict = {'keep_filtered': 0, 'max_gc': 60.0, 'min_gc': 40.0,
                     'min_length': 0, 'type_er_fl': 0, 'error_state': False,
                     'output_base_name': 'test', 'input_file': 'test.fastq',
                     'output_base_name_passed': 'test__passed.fastq',
                     'output_base_name_failed': 'test__failed.fastq'}
        # Correct gc bounds
        arg = ["fastq_updated_filter.py", "--gc_bounds", "40", "60", "test.fastq"]
        self.assertDictEqual(test_dict, arg_pars(arg))

        # Min bound > max. Changing positions
        arg = ["fastq_updated_filter.py", "--gc_bounds", "60", "40", "test.fastq"]
        self.assertDictEqual(test_dict, arg_pars(arg))

        # Type error (using default)
        arg = ["fastq_updated_filter.py", "--gc_bounds", "40w", "60b", "test.fastq"]
        test_dict["min_gc"] = 0.0
        test_dict["max_gc"] = 100.0
        test_dict["type_er_fl"] = 1
        self.assertDictEqual(test_dict, arg_pars(arg))

    def test_combined_args(self):
        test_dict = {'keep_filtered': 1, 'max_gc': 60.0, 'min_gc': 40.0,
                     'min_length': 20, 'type_er_fl': 0, 'error_state': False,
                     'output_base_name': 'test', 'input_file': 'test.fastq',
                     'output_base_name_passed': 'test__passed.fastq',
                     'output_base_name_failed': 'test__failed.fastq'}
        # Correct combined arguments
        arg = ["fastq_updated_filter.py", "--keep_filtered", "--min_length", "20", "--gc_bounds", "40", "60",
               "test.fastq"]
        self.assertDictEqual(test_dict, arg_pars(arg))

        # Some problem in args (no min_length gc bounds next)
        arg = ["fastq_updated_filter.py", "--keep_filtered", "--min_length", "--gc_bounds", "40", "60",
               "test.fastq"]
        test_dict["min_length"] = 0
        self.assertDictEqual(test_dict, arg_pars(arg))

    def test_arg_order(self):
        test_dict = {'keep_filtered': 1, 'max_gc': 60.0, 'min_gc': 40.0,
                     'min_length': 20, 'type_er_fl': 0, 'error_state': False,
                     'output_base_name': 'test', 'input_file': 'test.fastq',
                     'output_base_name_passed': 'test__passed.fastq',
                     'output_base_name_failed': 'test__failed.fastq'}
        diff_ordered_args = [
            ["fastq_updated_filter.py", "--keep_filtered", "--min_length", "20", "--gc_bounds", "40", "60",
             "test.fastq"],
            ["fastq_updated_filter.py", "--min_length", "20", "--keep_filtered", "--gc_bounds", "40", "60",
             "test.fastq"],
            ["fastq_updated_filter.py", "--min_length", "20", "--gc_bounds", "40", "60", "--keep_filtered",
             "test.fastq"],
            ["fastq_updated_filter.py", "--keep_filtered", "--gc_bounds", "40", "60", "--min_length", "20",
             "test.fastq"]
        ]
        for arg_set in diff_ordered_args:
            self.assertDictEqual(arg_pars(arg_set), test_dict)


class TestFastqFilter(unittest.TestCase):
    def setUp(self):
        self.default_test_arg = {'keep_filtered': 0, 'max_gc': 100.0, 'min_gc': 0.0,
                                 'min_length': 0, 'type_er_fl': 0, 'error_state': False,
                                 'output_base_name': 'test', 'input_file': 'test.fastq',
                                 'output_base_name_passed': 'test__passed.fastq',
                                 'output_base_name_failed': 'test__failed.fastq'}

    def test_out_file_exist(self):
        filter_reads(self.default_test_arg)
        self.assertTrue(os.path.exists("test__failed.fastq"))
        self.assertTrue(os.path.exists("test__passed.fastq"))
        os.remove("test__passed.fastq")
        os.remove("test__failed.fastq")

    def test_correct_gc_passed(self):
        # In test file there are 2 reads with 100 gc content
        new_gc_bound_arg = deepcopy(self.default_test_arg)
        new_gc_bound_arg["min_gc"] = 99
        filter_reads(new_gc_bound_arg)
        with open("test__passed.fastq") as passed_reads:
            length = 0
            names = []
            for line in passed_reads:
                length += 1
                if (length - 1) % 4 == 0:
                    names.append(line.strip())
        self.assertEqual(length / 4, 2)
        self.assertEqual(['@SRR13632571', '@SRR13632572'], names)

    def test_gc_both_bounds(self):
        # There is only one read with 50% gc
        new_gc_bound_arg = deepcopy(self.default_test_arg)
        new_gc_bound_arg["min_gc"] = 50
        new_gc_bound_arg["max_gc"] = 50
        filter_reads(new_gc_bound_arg)
        with open("test__passed.fastq") as passed_reads:
            length = 0
            names = []
            for line in passed_reads:
                length += 1
                if (length - 1) % 4 == 0:
                    names.append(line.strip())
        self.assertEqual(length / 4, 1)
        self.assertEqual(['@SRR13632573'], names)

    def test_random_gc_len(self):
        # There are 561 reads with such params
        new_gc_bound_arg = deepcopy(self.default_test_arg)
        new_gc_bound_arg["min_gc"] = 40
        new_gc_bound_arg["max_gc"] = 60
        new_gc_bound_arg["min_length"] = 98
        filter_reads(new_gc_bound_arg)
        with open("test__passed.fastq") as failed_file:
            s = failed_file.readlines()
        self.assertEqual(len(s)/4, 561)

    def test_no_failed_reads(self):
        # No filtered reads
        new_arg = deepcopy(self.default_test_arg)
        filter_reads(new_arg)
        with open("test__failed.fastq") as failed_file:
            s = failed_file.readlines()
        self.assertEqual(len(s), 0)
        with open("test__passed.fastq") as passed_file:
            s = passed_file.readlines()
        self.assertEqual(len(s) / 4, 678)

        # Save filtered reads
        new_arg = deepcopy(self.default_test_arg)
        new_arg["keep_filtered"] = 1
        new_arg["min_gc"] = 99
        filter_reads(new_arg)
        with open("test__failed.fastq") as failed_file:
            s = failed_file.readlines()
        self.assertGreater(len(s), 0)

    def test_length(self):
        # There are 3 reads shorter than 100
        new_arg = deepcopy(self.default_test_arg)
        new_arg["keep_filtered"] = 1
        new_arg["min_length"] = 100
        filter_reads(new_arg)
        with open("test__failed.fastq") as failed_file:
            s = failed_file.readlines()
        self.assertEqual(3, len(s)/4)

    def test_all_reads_in_out_files(self):
        # Check if there is no skiped reads (n_reads in input and output are equal)
        new_arg = deepcopy(self.default_test_arg)
        new_arg["keep_filtered"] = 1
        new_arg["min_length"] = 100
        new_arg["min_gc"] = 40
        new_arg["max_gc"] = 60
        filter_reads(new_arg)
        with open("test__passed.fastq") as passed_file:
            s = len(passed_file.readlines())
        with open("test__failed.fastq") as failed_file:
            s1 = len(failed_file.readlines())
        with open("test.fastq") as file:
            s_full = len(file.readlines())
        self.assertEqual((s+s1), s_full)

