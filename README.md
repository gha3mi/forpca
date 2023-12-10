[![GitHub](https://img.shields.io/badge/GitHub-ForPCA-blue.svg?style=social&logo=github)](https://github.com/gha3mi/forpca)
[![Version](https://img.shields.io/github/release/gha3mi/forpca.svg)](https://github.com/gha3mi/forpca/releases/latest)
[![Documentation](https://img.shields.io/badge/ford-Documentation%20-blueviolet.svg)](https://gha3mi.github.io/forpca/)
[![License](https://img.shields.io/github/license/gha3mi/forpca?color=green)](https://github.com/gha3mi/forpca/blob/main/LICENSE)
[![Build](https://github.com/gha3mi/forpca/actions/workflows/ci.yml/badge.svg)](https://github.com/gha3mi/forpca/actions/workflows/ci.yml)

<img alt="ForPCA" src="https://github.com/gha3mi/forpca/raw/main/media/logo.png" width="750">

**ForPCA**: A Fortran library for principal component analysis (PCA).

## Usage

```Fortran
use forpca, only: tpca
type(tpca) :: p

call p%pca(matrix, npc, method, coeff, score, latent, explained, matrix_app)

call p%finalize() ! finalize
```

## fpm dependency

If you want to use `ForPCA` as a dependency in your own fpm project,
you can easily include it by adding the following line to your `fpm.toml` file:

```toml
[dependencies]
forpca = {git="https://github.com/gha3mi/forpca.git"}
```

## How to run tests

**Reuirements:**

Fortran Compiler, LAPACK or MKL Libraries

**Clone the repository:**

You can clone the `ForPCA` repository from GitHub using the following command:

```shell
git clone https://github.com/gha3mi/forpca.git
```

```shell
cd forpca
```

**Intel Fortran Compiler (ifort)**

```shell
fpm @ifort-test
```

```shell
fpm @ifort-test-coarray
```

**Intel Fortran Compiler (ifx)**

```shell
fpm @ifx-test
```

```shell
fpm @ifx-test-coarray
```

**GNU Fortran Compiler (gfortran)**

```shell
fpm @gfortran-test
```

**NVIDIA Compiler (nvfortran)**

```shell
fpm @nvfortran-test
```

## API documentation

The most up-to-date API documentation for the master branch is available
[here](https://gha3mi.github.io/forpca/).
To generate the API documentation for `ForPCA` using
[ford](https://github.com/Fortran-FOSS-Programmers/ford) run the following
command:

```shell
ford ford.yml
```

## Contributing

Contributions to `ForPCA` are welcome!
If you find any issues or would like to suggest improvements, please open an issue.
