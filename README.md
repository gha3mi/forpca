![ForPCA](media/logo.png)
============

**ForPCA**: A Fortran library for principal component analysis (PCA).


## Installation

### fpm
ForPCA can be cloned and then built using [fpm](https://github.com/fortran-lang/fpm), following the instructions provided in the documentation available on Fortran Package Manager.

```bash
git clone https://github.com/gha3mi/forpca.git
cd forpca
fpm install --prefix .
```

Or you can easily include this package as a dependency in your `fpm.toml` file.

```toml
[dependencies]
forpca = {git="https://github.com/gha3mi/forpca.git"}
```
-----

## Documentation
To generate the documentation for the `ForPCA` module using [ford](https://github.com/Fortran-FOSS-Programmers/ford) run the following command:
```bash
ford ford.yml
```
-----

## Contributing
Contributions to `ForPCA` are welcome! If you find any issues or would like to suggest improvements, please open an issue or submit a pull request.