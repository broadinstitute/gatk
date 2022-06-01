import argparse


def scale_xy_bed_values(input_file, output_file, x_scale_factor, y_scale_factor):
    error = ""

    if x_scale_factor < 1.0:
        error = f"Error: illegal X chromosome weight scale factor {x_scale_factor}; scale factor must be >= 1.0\n"

    if y_scale_factor < 1.0:
        error = f"{error}Error: illegal Y chromosome weight scale factor {y_scale_factor}; scale factor must be >= 1.0"

    if len(error) > 0:
        raise ValueError(error)

    with open(input_file, 'r') as input_bed, open(output_file, 'w') as output_bed:
        while True:
            line = input_bed.readline()
            if not line:
                break

            line = line.rstrip('\n')

            if line.startswith('chrX') or line.startswith('chrY'):
                scale_factor = x_scale_factor if line.startswith('chrX') else y_scale_factor
                fields = line.split('\t')
                if len(fields) != 5:
                    raise ValueError(f"Expected 5 fields in input BED file, got {len(fields)}: {line}")
                weight = fields[-1]
                fields[-1] = str(int(int(weight) * scale_factor))
                output_bed.write('\t'.join(fields) + '\n')
            else:
                output_bed.write(line + '\n')


def parse_args():
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Scale X and Y BED values for more uniform extract shard runtimes')
    parser.add_argument('--input', type=str, help='Input BED file', required=True)
    parser.add_argument('--output', type=str, help='Output BED file', required=True)
    parser.add_argument('--xscale', type=float, help='X chromosome scaling factor', required=True)
    parser.add_argument('--yscale', type=float, help='Y chromosome scaling factor', required=True)

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    scale_xy_bed_values(args.input, args.output, args.xscale, args.yscale)
