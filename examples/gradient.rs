use std::{env,
          io::{BufWriter, Write},
          fs::File,
          error::Error};
use color_brewery::RGBColor;
use rgb::RGB8;
use color_brewery::{Palette, ColorRange};

type Err = Box<dyn Error>;

fn css_string(c: RGB8) -> String {
    format!("#{:02x}{:02x}{:02x}", c.r, c.g, c.b)
}

fn table_of_colors(fh: &mut impl Write, colors: &Vec<RGB8>,
                   width: u32, comment: &str) -> Result<(), Err> {
    writeln!(fh, "<table style=\"border: 0px;  border-spacing: 0px\"><tr>")?;
    for &c in colors {
        writeln!(fh, "  <td style=\"width: {width}px; height: 30px; \
                      background-color: {}\"></td>",
                css_string(c))?;
    }
    writeln!(fh, "<td rowspan=\"2\" style=\"padding-left: 7px\">\
                  {comment}</td></tr><tr>")?;
    for &c in colors {
        let c = c.to_gray();
        writeln!(fh, "  <td style=\"width: {width}px; height: 12px; \
                      background-color: {}\"></td>",
                 css_string(c))?;
    }
    writeln!(fh, "</tr></table><br/>")?;
    Ok(())
}

fn range(fh: &mut impl Write,
         color: impl Fn(f64) -> RGB8, n: usize,
         width: u32, comment: &str) -> Result<(), Err> {
    let dt = 1. / (n - 1) as f64;
    let colors: Vec<_> = (0 .. n).map(move |i| color(i as f64 * dt))
        .collect();
    table_of_colors(fh, &colors, width, comment)
}

fn gradient(fh: &mut impl Write, c0: [u8; 3], c1: [u8; 3], n: usize,
            width: u32, comment: &str) -> Result<(), Err> {
    let c0 = RGB8{ r: c0[0], g: c0[1], b: c0[2] };
    let c1 = RGB8{ r: c1[0], g: c1[1], b: c1[2] };
    let g = c0.gradient(&c1);
    range(fh, |t| g.rgb(t), n, width, comment)
}

fn palette(fh: &mut impl Write, palette: Palette<RGB8>, n: usize,
           width: u32, comment: &str) -> Result<(), Err> {
    let g = palette.gradient();
    range(fh, |t| g.rgb(t), n, width, comment)
}


fn main() -> Result<(), Err> {
    let mut fh = BufWriter::new(File::create("gradient.html")?);
    writeln!(fh, "<html>\n\
                  <head>\n\
                  <title>Color_brewery: test {}</title>\n\
                  </head>\n\
                  <body>",
             env::args().next().unwrap())?;
    writeln!(fh, "<h3>Hue</h3>")?;
    range(&mut fh, |t| RGB8::HUE.rgb(t), 10, 43, "")?;
    range(&mut fh, |t| RGB8::HUE.rgb(t), 30, 13, "")?;
    range(&mut fh, |t| RGB8::HUE.rgb(t), 150, 1, "")?;

    writeln!(fh, "<h3>Gradients</h3>")?;
    gradient(&mut fh, [94, 0, 99], [255, 235, 170], 10, 43, "")?;
    gradient(&mut fh, [94, 0, 99], [255, 235, 170], 30, 13, "")?;
    gradient(&mut fh, [94, 0, 99], [255, 235, 170], 150, 1, "")?;
    gradient(&mut fh, [255, 0, 0], [0, 0, 255], 150, 1, "")?;
    gradient(&mut fh, [255, 0, 0], [0, 255, 0], 150, 1,
             "Best <a href=\"https://youtu.be/XjHzLUnHeM0?t=230\"
              >to avoid red and green</a>.")?;
    gradient(&mut fh, [0, 0, 0], [255, 255, 255], 150, 1, "")?;
    gradient(&mut fh, [0, 0, 0], [0, 0, 255], 150, 1, "")?;
    gradient(&mut fh, [0, 0, 128], [144, 144, 255], 150, 1, "")?;

    writeln!(fh, "<h3>Palettes (sequential and diverging)</h3>")?;
    for (comment, p) in [("viridis", RGB8::viridis()),
                         ("magma", RGB8::magma()),
                         ("inferno", RGB8::inferno()),
                         ("plasma", RGB8::plasma())] {
        palette(&mut fh, p, 256, 1, comment)?;
    }

    macro_rules! palette { ($($name: ident),*) => {
        $( if let Some(p) = RGB8::$name().last() {
            let c = format!("{} ({} colors)", stringify!($name), p.len());
            table_of_colors(&mut fh, &p.colors(), 40, &c)?;
            palette(&mut fh, p, 128, 1,
                    &format!("{} (interpolated)", stringify!($name)))?;
        } )*
    }}
    palette!(ylgn, ylgnbu, gnbu, bugn, pubugn, pubu, bupu, rdpu, purd,
             orrd, ylorrd, ylorbr, purples, blues, greens, oranges,
             reds, greys,
             // Diverging
             puor, brbg, prgn, piyg, rdbu, rdgy, rdylbu, spectral, rdylgn);

    writeln!(fh, "<h3>Palettes (others)</h3>")?;
    macro_rules! palette { ($($name: ident),*) => {
        $( if let Some(p) = RGB8::$name().last() {
            let c = format!("{} ({} colors)", stringify!($name), p.len());
            table_of_colors(&mut fh, &p.colors(), 40, &c)?;
        } )*
    }}
    palette!(set1, pastel1, set2, pastel2, dark2, set3, paired, accent);

    writeln!(fh, "</body>\n\
                  </html>")?;
    Ok(())
}
