from . import arguments, train, io, plot


def main():
    args = arguments.parse_args()
    output = train.run(args)
    io.write_output(input_vcf_path=args.vcf, output_vcf_path=args.output, output_data=output)
    if args.plots_dir is not None:
        plot.plot_sites(vcf_path=args.output, out_dir=args.plots_dir)


if __name__ == '__main__':
    main()
