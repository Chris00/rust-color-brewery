// To use this script, download the JSON files
// - http://colorbrewer2.org/colorbrewer_schemes.js and
// - http://colorbrewer2.org/cmyk.js
// and copy them in tools/.  Then remove everything at be beginning
// before the first "{" and replace all single quotes by double quotes.

use std::{string::String,
          collections::BTreeMap,
          env,
          fs::File,
          io::{BufWriter, prelude::*},
          error::Error};
use serde_hjson::Value::{self, *};

mod matplotlib;

#[derive(Debug, Clone, Copy)]
enum Trivalent {
    Yes,
    Maybe,
    No,
}

#[derive(Debug)]
struct Palette {
    rgb: Vec<[f64; 3]>,// RGB Colors of the palette
    typ: String,
    blind: Trivalent,
    print: Trivalent,
    photocopy: Trivalent,
    lcd: Trivalent,
}

type Palettes = BTreeMap<String, Vec<Palette>>;

fn capitalize(s: &str) -> String {
    s.chars().enumerate().map(|(i, c)| {
        if i == 0 { c.to_ascii_uppercase() }
        else { c }}).collect()
}

/// Convert colors such as "rgb(67,147,195)" to `[67., 147., 195.]`.
fn parse_color(s: &str) -> [f64; 3] {
    let mut colors = [0.; 3];
    for (i, c) in s[4 .. s.len() - 1].split(',').enumerate() {
        let c = match str::parse::<f64>(c) {
            Ok(c) => c,
            Err(_) => panic!("color-brewery-tools: not a color “{}”", s) };
        colors[i] = c
    }
    colors
}

fn trivalent_of_int(i: &Value) -> Trivalent {
    match i {
        I64(0) | U64(0) => Trivalent::No,
        I64(1) | U64(1) => Trivalent::Yes,
        I64(2) | U64(2) => Trivalent::Maybe,
        _ => panic!("color-brewery-tools: {i:?} not an int in 0..=2"),
    }
}

fn trivalent_prop(prop: Option<&Value>) -> Vec<Trivalent> {
    match prop {
        Some(Array(i)) => {
            if i.len() == 1 {
                vec![trivalent_of_int(&i[0]); 12]
            } else {
                let mut prop: Vec<_> = i.iter().map(trivalent_of_int).collect();
                // Some property lists are not long enough, assume
                // when not set, the last item is always "No".
                prop.push(Trivalent::No);
                prop
            }
        }
        None => vec![Trivalent::Maybe; 12],
        _ => panic!("color-brewery-tools: {prop:?}"),
    }
}

fn add_palettes_from_json(map: &mut Palettes, json: Value) {
    if let Object(m) = json {
        for (name, palettes) in m.into_iter() {
            let map = map.entry(name).or_default();
            if let Object(palettes) = palettes {
                if let Some(Object(p)) = palettes.get("properties") {
                    let ty = match p.get("type") {
                        Some(String(t)) => t,
                        _ => panic!("color-brewery-tools") };
                    let blind = trivalent_prop(p.get("blind"));
                    let print = trivalent_prop(p.get("print"));
                    let photocopy = trivalent_prop(p.get("copy"));
                    let lcd = trivalent_prop(p.get("screen"));
                    let rgbs = palettes.iter()
                        .filter(|(n,_)| n.as_str() != "properties");
                    for (i, (n, rgb)) in rgbs.enumerate() {
                        let n = str::parse::<usize>(n).unwrap();
                        assert_eq!(i + 3, n);
                        let rgb = match rgb {
                            Array(a) => a,
                            _ => panic!("color-brewery-tools") };
                        let rgb: Vec<_> = rgb.iter()
                            .map(|v| match v {
                                String(c) => parse_color(c),
                                _ => panic!("color-brewery-tools") })
                            .collect();
                        let palette = Palette {
                            rgb,
                            typ: ty.to_string(),
                            blind: blind[i],
                            print: print[i],
                            photocopy: photocopy[i],
                            lcd: lcd[i],
                        };
                        map.push(palette)
                    }
                }
            }
        }
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut brewer_maps = Palettes::new();

    let fh_rgb = File::open("colorbrewer_schemes.js")?;
    let rgb_json: Value = serde_hjson::from_reader(fh_rgb)?;
    add_palettes_from_json(&mut brewer_maps, rgb_json);

    let mut fh = BufWriter::new(File::create("../src/palettes.rs")?);
    writeln!(fh, "// Written by {}\n\n\
                  use lazy_static::lazy_static;\n\
                  pub(crate) mod ty;\n\
                  use ty::*;",
             env::args().next().unwrap())?;

    writeln!(fh, "\n// Matplotlib palettes")?;
    fn matplotlib<const N: usize>(
        fh: &mut impl Write,
        name: &str,
        palette: [[f64; 3]; N]) -> std::io::Result<()> {
        write!(fh, "lazy_static! {{\n  \
                    pub(crate) static ref {name}: PaletteData = {{\n  \
                    PaletteData {{\n    \
                    typ: PaletteType::Seq,\n    \
                    blind: Trivalent::Yes,\n    \
                    print: Trivalent::Yes,\n    \
                    photocopy: Trivalent::Maybe,\n    \
                    lcd: Trivalent::Yes,\n    \
                    rgb: vec![\n")?;
        for [r, g, b] in palette {
            writeln!(fh, "      RGBA{{r: {:>10.6}, g: {:>10.6}, b: {:>10.6}, \
                          a: 1.}},", 255. * r, 255. * g, 255. * b)?;
        }
        writeln!(fh, "    ]}}\n  }};\n}}\n")?;
        Ok(())
    }
    matplotlib(&mut fh, "MAGMA", matplotlib::MAGMA_RGB)?;
    matplotlib(&mut fh, "INFERNO", matplotlib::INFERNO_RGB)?;
    matplotlib(&mut fh, "PLASMA", matplotlib::PLASMA_RGB)?;
    matplotlib(&mut fh, "VIRIDIS", matplotlib::VIRIDIS_RGB)?;
    matplotlib(&mut fh, "CIVIDIS", matplotlib::CIVIDIS_RGB)?;
    matplotlib(&mut fh, "TWILIGHT", matplotlib::TWILIGHT_RGB)?;
    matplotlib(&mut fh, "TURBO", matplotlib::TURBO_RGB)?;

    writeln!(fh, "// Brewer colormaps — see http://colorbrewer2.org/\n\
                  // Number of maps: {}",
             brewer_maps.len())?;
    for (name, palettes) in &brewer_maps {
        writeln!(fh, "lazy_static! {{\n  \
                      pub(crate) static ref {}: Vec<PaletteData> = {{\n    \
                      vec![",
                 name.to_ascii_uppercase())?;
        for p in palettes {
            write!(fh, "      PaletteData {{\n        \
                        typ: PaletteType::{},\n        \
                        blind: Trivalent::{:?},\n        \
                        print: Trivalent::{:?},\n        \
                        photocopy: Trivalent::{:?},\n        \
                        lcd: Trivalent::{:?},\n        \
                        rgb: vec![",
                   capitalize(&p.typ), p.blind, p.print, p.photocopy, p.lcd)?;
            for [r, g, b] in &p.rgb {
                write!(fh, "RGBA{{r: {r:.1}, g: {g:.1}, b: {b:.1}, a: 1.}},")?;
            }
            writeln!(fh, "]\n      }},")?;
        }
        writeln!(fh, "      ]\n  }};\n}}\n")?;
    }

    write!(fh, "lazy_static! {{\n  \
                pub(crate) static ref ALL_MATPLOTLIB_PALETTES: \
                [&'static PaletteData; 4] = [\n    \
                MAGMA.deref(), INFERNO.deref(), PLASMA.deref(), \
                VIRIDIS.deref()];\n  \
                pub(crate) static ref ALL_BREWER_PALETTES: \
                [&'static Vec<PaletteData>; {}] = {{\n    [",
           brewer_maps.len())?;
    for (name, _) in &brewer_maps {
        write!(fh, "{}.deref(),\n     ", name.to_ascii_uppercase())?;
    }
    writeln!(fh, "]\n  }};\n}}")?;

    Ok(())
}
