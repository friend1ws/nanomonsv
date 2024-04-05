#! /usr/bin/env python3

import unittest
import subprocess

class TestMain(unittest.TestCase):

    def setUp(self):
        pass

    def test_command(self):
        subprocess.run(["nanomonsv", "--version"])
        

if __name__ == "__main__":
    unittest.main()

