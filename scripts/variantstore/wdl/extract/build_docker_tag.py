from datetime import datetime
from subprocess import run
import argparse
import subprocess


def build_tag(args):
    if not args.release and not args.branch:
        raise ValueError("Neither release nor branch option specified.")
    if args.release and args.branch:
        raise ValueError("Both release and branch options specified.")

    current_date = datetime.now()
    date_string = current_date.strftime('%Y-%m-%d')

    if args.release:
        return f"{date_string}-alpine"
    else:
        if not args.dummy_testing_hash:
            proc = run(["git", "rev-parse", "--short", "HEAD"], stdout=subprocess.PIPE)
            git_hash = proc.stdout.rstrip().decode('utf-8')
        else:
            # Do not actually try to run git during testing, the .git directory is not mounted into the container.
            git_hash = args.dummy_testing_hash
        return f"{date_string}-alpine-{git_hash}"


def build_argument_parser():
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Build a tag for Variants Docker image')
    parser.add_argument('-r', '--release', action='store_true', help='Create a release tag')
    parser.add_argument('-b', '--branch', action='store_true', help='Create a branch tag')
    parser.add_argument('-d', '--dummy-testing-hash', default=None,
                        help='Dummy short git hash to return during testing')
    return parser


if __name__ == '__main__':
    arg_parser = build_argument_parser()
    parsed_args = arg_parser.parse_args()
    tag = build_tag(parsed_args)
    print(tag, end='')
