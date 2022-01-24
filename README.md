# prot-scriber

Assigns short human readable descriptions (HRD) to query biological sequences using reference canditate descriptions. In this, `prot-scriber` consumes sequence similarity search (Blast or Diamond or similar) results in tabular format. A customized lexical analysis is carried out on the descriptions of these Blast Hits and a resulting HRD is assigned to the query sequences. 

`prot-scriber` can also apply the same methodology to produce HRDs for sets of biological sequences, i.e. gene families. 

## Installation

### Pre built binaries

You can choose to download a pre-built binary ready to be executed on your platform and operating system (see table below) from the [Releases page](https://github.com/usadellab/prot-scriber/releases). We strongly recommend to _rename_ the downloaded release file to a simple `prot-scriber` (or `prot-scriber.exe` on Windows).

|Operating System|CPU-Architecture|Release-Name|Comment|
|---|---|---|---|
|Windows 7 or higher|any x86|windows_prot-scriber.exe|to be used in a terminal (`cmd` or Power-Shell)|
|any GNU-Linux|any Intel x86, 64 bits|x86_64-unknown-linux-gnu_prot-scriber||
|any GNU-Linux|any aarch, 64 bits|aarch64-unknown-linux-gnu_prot-scriber|e.g. for Raspberry Pi|
|Apple / Mac OS|any Mac Computer with Mac OS 10|x86_64-apple-darwin_prot-scriber||

### Compilation from source code

#### Prerequisites 

`prot-scriber` is written in Rust. That makes it extremely performant and, once compiled for your operating system (OS), can be used on any machine with that particular OS environment. To compile it for your platform you first need have Rust and cargo installed. Follow [the official instructions](https://www.rust-lang.org/tools/install).

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

The above compilation command has generated an executable binary file in the current directoy `./target/release/prot-scriber`. You can just go ahead and use `prot-scriber` now that you have compiled it successfully (see Usage below).

#### Global installation

If you are familiar with installing self compiled tools on a system wide level, this section will provide no news to you. It is convenient to make the compiled executable `prot-scriber` programm available from anywhere on your system. To achieve this, you need to copy it to any place you typically have your programms installed, or add its directory to your, our all users' `$PATH` environment. In doing so, e.g. in case you are a system administrator, you make `prot-scriber` available for all users of your infrastructe. Make sure you and your group have executable access rights to the file. You can adjust these access right with `chmod ug+x ./target/release/prot-scriber`. You, and possibly other users of your system, are now ready to run `prot-scriber`.

## Usage

`prot-scriber` is a command line tool and _must_ be used in a terminal application. On Windows that will be `cmd` or PowerShell, on Mac OS X or any Linux / Unix system that will be a standard terminal shell.

To use `prot-scriber` please read the manual. All your questions will be answered there. Please use 
```sh
prot-scriber --help
```
to read it.

_Happy `prot-scribing`!_
