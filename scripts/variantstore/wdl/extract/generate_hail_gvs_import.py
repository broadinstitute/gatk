import json
from collections import defaultdict


def generate_avro_dict(avro_prefix, listing):
    avro_file_arguments = defaultdict(list)
    super_partitioned_keys = {'vets', 'refs'}
    slashed_avro_prefix = avro_prefix + "/" if not avro_prefix[-1] == '/' else avro_prefix

    with open(listing, 'r') as f:
        for full_path in f.readlines():
            full_path = full_path.strip()
            if not full_path.endswith(".avro"):
                continue
            relative_path = full_path[len(slashed_avro_prefix):]
            parts = relative_path.split('/')
            key = parts[0]

            if key in super_partitioned_keys:
                # Get the zero based index from the `vet_001` component
                index = int(parts[1].split('_')[-1]) - 1
                if len(avro_file_arguments[key]) == index:
                    avro_file_arguments[key].append([])
                avro_file_arguments[key][index].append(full_path)
            else:
                avro_file_arguments[key].append(full_path)

    return dict(avro_file_arguments)


def generate_avro_args(avro_dict):
    args = []
    for k, v in avro_dict.items():
        s = f'{k}={json.dumps(v, indent=4)}'
        args.append(s)
    return ',\n'.join(args)
