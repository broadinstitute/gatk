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
        proc = run(["git", "rev-parse", "--short", "HEAD"], stdout=subprocess.PIPE)
        git_hash = proc.stdout.rstrip().decode('utf-8')
        return f"{date_string}-alpine-{git_hash}"


def build_argument_parser():
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Build a tag for Variants Docker image')
    parser.add_argument('-r', '--release', action='store_true', help='Create a release tag')
    parser.add_argument('-b', '--branch', action='store_true', help='Create a branch tag')
    return parser


if __name__ == '__main__':
    arg_parser = build_argument_parser()
    tag = build_tag(arg_parser.parse_args())
    print(tag, end='')
