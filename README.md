Color brewery
=============

This crate define colors palettes an provide a searchable access to
find those matching some criteria.  In particular, it contains all the
palettes of Matplotlib and of [ColorBrewer][] designed by the
cartographer Cynthia A. Brewer to assist mapmakers in choosing
appropriate color schemes.  As a convenience, color gradients based on
those palettes are also implemented.

The default representation of colors is the one given by the [rgb
crate][] but any type can be used: it suffices it implements the
`RGBColor` trait.

[ColorBrewer]: http://colorbrewer2.org/
[rgb crate]: https://crates.io/crates/rgb

## Installation

Run the following Cargo command in your project directory:
```
cargo add color-brewery
```

or add the following line to your `Cargo.toml`:

```
[dependencies]
color-brewery = "0.1"
```

## Documentation

See [doc.rs](https://docs.rs/color-brewery).
