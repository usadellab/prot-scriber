# prot-scriber

Assigns short human readable descriptions (HRD) to query biological sequences using reference candidate descriptions. In this, `prot-scriber` consumes sequence similarity search (Blast or Diamond or similar) results in tabular format. A customized lexical analysis is carried out on the descriptions of these Blast Hits and a resulting HRD is assigned to the query sequences. 

`prot-scriber` can also apply the same methodology to produce HRDs for sets of biological sequences, i.e. gene families. 

## Installation

### Pre built binaries

You can choose to download a pre-built binary ready to be executed from the table below, if you want the latest stable version. Other versions can be downloaded from the [Releases page](https://github.com/usadellab/prot-scriber/releases). Have a look at the below table to know which is the version you need for your operating system and platform. 

We strongly recommend to _rename_ the downloaded release file to a simple `prot-scriber` (or `prot-scriber.exe` on Windows).

Note that on Mac OS and Unix / Linux operating systems you need to make the downloaded and renamed `prot-scriber` file _executable_. In order to achieve this, open a terminal shell, navigate (`cd`) to the directory where you saved `prot-scriber`, and execute `chmod a+x ./prot-scriber`.


|Operating System|CPU-Architecture|Release-Name (click to download)|Comment|
|---|---|---|---|
|Windows 7 or higher|any|[windows_prot-scriber.exe](https://github.com/usadellab/prot-scriber/releases/download/latest-stable/x86_64-pc-windows-gnu_prot-scriber.exe)|to be used in a terminal (`cmd` or Power-Shell)|
|any GNU-Linux|any Intel x86, 64 bits|[x86_64-unknown-linux-gnu_prot-scriber](https://github.com/usadellab/prot-scriber/releases/download/latest-stable/x86_64-unknown-linux-gnu_prot-scriber)||
|any GNU-Linux|any aarch, 64 bits|[aarch64-unknown-linux-gnu_prot-scriber](https://github.com/usadellab/prot-scriber/releases/download/latest-stable/aarch64-unknown-linux-gnu_prot-scriber)|e.g. for Raspberry Pi|
|Apple / Mac OS|any Mac Computer with Mac OS 10|[x86_64-apple-darwin_prot-scriber](https://github.com/usadellab/prot-scriber/releases/download/latest-stable/x86_64-apple-darwin_prot-scriber)||

### Compilation from source code

#### Prerequisites 

`prot-scriber` is written in Rust. That makes it extremely performant and, once compiled for your operating system (OS), can be used on any machine with that particular OS environment. To compile it for your platform you first need to have Rust and cargo installed. Follow [the official instructions](https://www.rust-lang.org/tools/install).

#### Obtain the code

Download the [latest stable release of `prot-scriber` here](https://github.com/usadellab/prot-scriber/archive/refs/tags/latest-stable.zip).

Unzip it, e.g. by double clicking it or by using the command line:
```sh
unzip latest-stable.zip
```

#### Compile `prot-scriber`

Change into the directory of the downloaded `prot-scriber` code, e.g. in a Mac OS or Linux terminal `cd prot-scriber-latest-stable` after having unpacked the latest stable release (see above).

Now, compile `prot-scriber` with
```sh
cargo build --release
```

The above compilation command has generated an executable binary file in the current directory `./target/release/prot-scriber`. You can just go ahead and use `prot-scriber` now that you have compiled it successfully (see Usage below).

#### Global installation

If you are familiar with installing self compiled tools on a system wide level, this section will provide no news to you. It is convenient to make the compiled executable `prot-scriber` program available from anywhere on your system. To achieve this, you need to copy it to any place you typically have your programs installed, or add its directory to your, our all users' `$PATH` environment. In doing so, e.g. in case you are a system administrator, you make `prot-scriber` available for all users of your infrastructure. Make sure you and your group have executable access rights to the file. You can adjust these access right with `chmod ug+x ./target/release/prot-scriber`. You, and possibly other users of your system, are now ready to run `prot-scriber`.

## Usage

`prot-scriber` is a command line tool and _must_ be used in a terminal application. On Windows that will be `cmd` or PowerShell, on Mac OS X or any Linux / Unix system that will be a standard terminal shell.

To use `prot-scriber` please read the [**manual** (click here for the latest stable version)](https://github.com/usadellab/prot-scriber/releases/download/latest-stable/prot-scriber-manual.txt). All your questions will be answered there. 

In the command prompt (`cmd` or PowerShell on Windows) or the Terminal (Mac OS, Linux, or Unix) use  
```sh
prot-scriber --help
```
to get the manual.

_Happy `prot-scribing`!_

## Development / Contribute

`prot-scriber` is open source. Please feel free and invited to contribute. 

### Preparation of releases (pre-compiled executables)

This repository is set up to use [GitHub Actions](https://github.com/features/actions) (see `.github/workflows/push.yml` for details). We use GitHub actions to trigger compilation of `prot-scriber` every time a Git Tag is pushed to this repository. So, if you, after writing new code and committing it, do the following in your local repo:

```sh
git tag -a 'version-Foo-Bar-Baz' -m "My fancy new version called Foo Bar Baz"
git push origin master --tags
```

GitHub will automatically compile `prot-scriber` and provide executable binary versions of prot-scriber for the platforms and operating systems mentioned above. The resulting binaries are then made available for download on the [releases page](https://github.com/usadellab/prot-scriber/releases).

In short, you do not need to worry about compiling your latest version to make it available for download and the different platforms and operating systems, GitHub Actions take care of this for you. 
