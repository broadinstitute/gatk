import re
import unittest
from build_docker_tag import build_tag, build_argument_parser


class TestBuildDockerTag(unittest.TestCase):
    # noinspection PyPep8Naming
    def setUp(self) -> None:
        self.argument_parser = build_argument_parser()

    def test_tag(self):
        tag_re = re.compile(r"^20\d{2}-\d{2}-\d{2}-alpine-f00ba4ba5$")

        args = self.argument_parser.parse_args(['--dummy-testing-hash', 'f00ba4ba5'])
        tag = build_tag(args)
        self.assertTrue(tag_re.match(tag))
