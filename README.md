[![crates.io](https://img.shields.io/crates/v/ent-rs.svg)](https://crates.io/crates/ent-rs)
[![docs.rs](https://docs.rs/ent-rs/badge.svg)](https://docs.rs/ent-rs)
[![CI](https://github.com/chbtoys/ent-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/chbtoys/ent-rs/actions)

# ent-rs

**ent-rs** is a Rust library for analyzing the entropy and randomness characteristics of binary data. It is a Rust port and extension of the functionality from Tucker's `ent.hpp`.

It is useful for evaluating:
- Compression potential
- Statistical randomness
- Distribution uniformity
- Pseudorandom number generators
- Data quality or obfuscation

## Features

- Byte and bit mode entropy analysis
- Chi-square test and p-value
- Arithmetic mean
- Monte Carlo Pi estimation
- Serial correlation
- Value frequency tables

## Usage

Add to your `Cargo.toml`:

```toml
ent-rs = "0.1"
```

## Example

```rust
use ent_rs::EntStats;
let data = std::fs::read("mydata.bin").unwrap();
let stats = EntStats::from_data(&data, false);
println!("Entropy: {:.4}", stats.entropy);
```
