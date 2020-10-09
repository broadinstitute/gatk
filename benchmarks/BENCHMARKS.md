# Data generation benchmarks :straight_ruler:
Benchmarks keep track of how quickly we can produce data for training models.
To run all benchmarks:
```bash
python benchmarks/benchmark.py
```


## Index
* [Running benchmarks](#running-benchmarks)
* [Contributing benchmarks](#running-benchmarks)


## Running benchmarks
Benchmarks are run using [benchmarks.py](./benchmark.py).
You can specify specific benchmarks using
```bash
python benchmarks/benchmarks.py [name of benchmark] [name of another benchmark] ...
```
The names of the currently available benchmarks are the directory names under [benchmark_results](./benchmark_results).

## Contributing benchmarks

There are three components to a benchmark:
1. [The type of data](#data-descriptions)
2. [The way the data is iterated over](#generatorfactories)

All three can be tinkered with.

### Data descriptions
Synthetic data is produced in [data.py](./data.py) using the function `data.build_example`.
A synthetic datum is described by a 3-tuple
```
[name, shape, data type]
```
For example, the `ecg_single_task` benchmark simulates reading ECG bmi pairs,
so it has two data descriptions:
```python
('ecg', (5000, 12), StorageType.CONTINUOUS),
('bmi', (1,), StorageType.CONTINUOUS),
```

### GeneratorFactories
`GeneratorFactorie`s prepare the data for a class of generators and produce data generators that use that data.
for example `benchmark.TensorGeneratorFactory` builds `ml4h.TensorGenerator`s with `hd5`s to store data.
`GeneratorFactorie`s must implement a `setup` function which prepares the synthetic data given a `DataDescription`
and a number of synthetic samples to build.

`GeneratorFactorie`s must also implement `__call__(batch_size: int, num_workers: int)`
which produces a data generator with that batch size and number of multiprocessing workers.
