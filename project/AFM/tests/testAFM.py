from unittest import TestCase

import AFM
from AFM.command_line import main

class TestJoke(TestCase):
    def test_is_string(self):
        s = AFM.joke()
        self.assertTrue(isinstance(s, basestring))

class TestConsole(TestCase):
    def test_basic(self):
        main()
