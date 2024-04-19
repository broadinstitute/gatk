from datetime import datetime
from subprocess import run
import argparse
import subprocess


def build_tag(args):
    current_date = datetime.now()
    date_string = current_date.strftime('%Y-%m-%d')

    if not args.dummy_testing_image_id:
        proc = run(["bash", "-c", "docker images --quiet | head -1"], stdout=subprocess.PIPE)
        docker_image_id = proc.stdout.rstrip().decode('utf-8')
    else:
        # Do not actually try to run git during testing, the .git directory is not mounted into the container.
        docker_image_id = args.dummy_testing_image_id
    return f"{date_string}-alpine-{docker_image_id}"


def build_argument_parser():
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Build a tag for Variants Docker image')
    parser.add_argument('-d', '--dummy-testing-image-id', default=None,
                        help='Dummy short Docker image ID to return during testing')
    return parser


if __name__ == '__main__':
    arg_parser = build_argument_parser()
    parsed_args = arg_parser.parse_args()
    tag = build_tag(parsed_args)
    print(tag, end='')
