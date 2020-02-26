from . import train, arguments, plot, infer, io


def train():
    args = arguments.parse_args()
    train.run(args)


def infer():
    args = arguments.parse_args()
    output, global_stats = infer.run(args)
    io.write_vcf(input_vcf_path=args.vcf, output_vcf_path=args.output, output_data=output, global_stats=global_stats)
    if args.plots_dir is not None:
        plot.plot_sites(vcf_path=args.output, out_dir=args.plots_dir)


def main():
    train()


if __name__ == '__main__':
    main()
