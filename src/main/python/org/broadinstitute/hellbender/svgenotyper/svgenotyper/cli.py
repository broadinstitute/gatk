from . import arguments, train


def main():
    args = arguments.parse_args()
    print(train.run(args))


if __name__ == '__main__':
    main()
