from . import arguments, train, io


def main():
    args = arguments.parse_args()
    output = train.run(args)
    io.write_output(input_vcf_path=args.vcf, output_vcf_path=args.output, output_data=output)


if __name__ == '__main__':
    main()
