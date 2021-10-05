import argparse

def generate_ordered_paths(root_path, path_suffix, number):
    for i in range(0, number - 1):
        print(f"{root_path}{i}{path_suffix}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Extract subpopulation per sample data out of a callset TSV')
    parser.add_argument('--root_path',type=str, metavar='string', help='path plus file prefix', required=True)
    parser.add_argument('--path_suffix',type=str, metavar='integer', help='path suffix', required=True)
    parser.add_argument('--number',type=int, metavar='string', help='number of files', required=True)

    args = parser.parse_args()

    generate_ordered_paths(args.root_path,
                          args.path_suffix,
                          args.number)
