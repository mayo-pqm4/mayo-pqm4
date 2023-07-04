This repository contains the code accompanying the paper **MAYO: Optimized Implementation with Revised Parameters for ARMv7-M**. The paper is available [here](https://eprint.iacr.org/2023/540).

Authors:
- Arianna Gringiani
- Alessio Meneghetti
- Edoardo Signorini
- Ruggero Susella
- Giovanni Tognolini

**Note**: This repository does not contain the official implementation of MAYO and is not related to MAYO's submission to NIST.

# Structure

We use the [pqm4](https://github.com/mupq/pqm4) framework for benchmarks and board integration.

The main C and ASM implementation of MAYO can be found in [`MAYO/`](./MAYO/). Relevant files are symlinked in the implementation folders in [`crypto_sign/`](./crypto_sign/).

# Results

We target the [NUCLEO-L4R5ZI](https://www.st.com/en/evaluation-tools/nucleo-l4r5zi.html) board with a STM32L4R5ZIT6 Cortex-M4 core, featuring 2MB of Flash and 640KB of RAM.

| Parameters     | Sign (K cycles) | Verify (K cycles) |
|----------------|-----------------|-------------------|
| mayo-orig-1    | 157,134         | 38,880            |
| mayo-new-1     | 50,183          | 7,371             |
| mayo-new-1-k14 | 126,849         | 35,700            |

In the paper we also target the STM32H753 Nucleo-144 board with a STM32H753ZIT6 Cortex-M7 core. This board is not supported by pqm4. The results for this board are available in Table 4 of the paper.

# Setup

Clone the repository and its submodules with the `--recursive` option:

```
git clone --recurse-submodule https://github.com/mayo-pqm4/mayo-pqm4
```

If already cloned, use:

```
git submodule update --init --recursive
```

Follow pqm4 [instructions](https://github.com/mupq/pqm4#setupinstallation) for framework and board setup. Verify proper functioning on the QEMU simulator of mps2-an386 board by testing the implementation of dilithium in the original [`pqm4/`](./pqm4/) folder:

```
cd pqm4
python3 test.py -p mps2-an386 -u /dev/ttyACM0 dilithium2
```

and with the nucleo-l4r5zi board:

```
python3 test.py -p nucleo-l4r5zi -u /dev/ttyACM0 dilithium2
```

# Benchmarking

We use the [`benchmarks.py`](./benchmarks.py) script from pqm4 for benchmarking. The parameters to be tested can be specified at the end of the command or left blank to test all parameter sets.
The code in this repository is optimized for signing and verification procedures only. We use a fork of pqm4 where the mupq speed benchmark performs key generation only once at the beginning of the test.

 The following options can be used:

- `-i`: selects the number of benchmark iterations
- `--nostack`: skip stack consumption benchmark
- `--nospeed`: skip speed benchmark
- `--nohashing`: skip the measure of cycles spent in symmetric primitives (SHA-2, SHA-3, and AES) evaluations
- `--nosize`: skip code size measurement

To test all parameter sets on nucleo-l4r5zi, use:

```
python3 benchmarks.py -p nucleo-l4r5zi -u /dev/ttyACM0 -o speed -i 10000
```

The benchmark results are stored in `benchmarks/`. They can be automatically converted to a markdown table or to csv using the [`convert_benchmarks.py`](./convert_benchmarks.py) script:

```
# to markdown
python3 convert_benchmarks.py md 

# to csv
python3 convert_benchmarks.py csv
```

# License

The implementation in [`MAYO/`](./MAYO/) is released under [CC0](https://creativecommons.org/publicdomain/zero/1.0/).
For third-party code we refer to the relevant licenses, indicated in the respective repositories.