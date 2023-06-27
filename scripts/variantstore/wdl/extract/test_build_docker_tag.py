import re
import unittest
from build_docker_tag import build_tag, build_argument_parser


class TestBuildDockerTag(unittest.TestCase):
    def setUp(self) -> None:
        self.argument_parser = build_argument_parser()

    def test_release(self):
        tag_re = re.compile(r"^20\d{2}-\d{2}-\d{2}-alpine$")

        args = self.argument_parser.parse_args(['-r'])
        self.assertTrue(args.release)
        self.assertFalse(args.branch)

        tag = build_tag(args)
        self.assertTrue(tag_re.match(tag))

        args = self.argument_parser.parse_args(['--release'])
        self.assertTrue(args.release)
        self.assertFalse(args.branch)

        tag = build_tag(args)
        self.assertTrue(tag_re.match(tag))

    def test_branch(self):
        tag_re = re.compile(r"^20\d{2}-\d{2}-\d{2}-alpine-[0-9a-f]{9}$")

        args = self.argument_parser.parse_args(['-b'])
        self.assertFalse(args.release)
        self.assertTrue(args.branch)

        tag = build_tag(args)
        self.assertTrue(tag_re.match(tag))

        args = self.argument_parser.parse_args(['--branch'])
        self.assertFalse(args.release)
        self.assertTrue(args.branch)

        tag = build_tag(args)
        self.assertTrue(tag_re.match(tag))

    def test_no_args(self):
        with self.assertRaises(ValueError) as context:
            args = self.argument_parser.parse_args([])
            build_tag(args)
        self.assertTrue('Neither release nor branch option specified.' in context.exception.args)

    def test_too_many_args(self):
        with self.assertRaises(ValueError) as context:
            args = self.argument_parser.parse_args(['--branch', '--release'])
            build_tag(args)
        self.assertTrue('Both release and branch options specified.' in context.exception.args)
