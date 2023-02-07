use std::{env,
          io::{BufWriter, Write},
          fs::File,
          ops::Deref,
          error::Error};
use colorsys::Rgb;
use color_brewery::{Hue, Gradient, ColorRange, Palette,
                    VIRIDIS, MAGMA, INFERNO, PLASMA};

type Err = Box<dyn Error>;

fn to_gray(c: &Rgb) -> Rgb {
  let x = 0.299 * c.red() + 0.587 * c.green() + 0.114 * c.blue();
    Rgb::from((x, x, x))
}

fn table_of_colors(fh: &mut impl Write, colors: &Vec<Rgb>,
                   width: u32, comment: &str) -> Result<(), Err> {
    writeln!(fh, "<table style=\"border: 0px;  border-spacing: 0px\"><tr>")?;
    for c in colors {
        writeln!(fh, "  <td style=\"width: {width}px; height: 30px; \
                      background-color: {}\"></td>",
                 c.to_css_string())?;
    }
    writeln!(fh, "<td rowspan=\"2\" style=\"padding-left: 7px\">\
                  {comment}</td></tr><tr>")?;
    for c in colors {
        let c = to_gray(c);
        writeln!(fh, "  <td style=\"width: {width}px; height: 12px; \
                      background-color: {}\"></td>",
                 c.to_css_string())?;
    }
    writeln!(fh, "</tr></table><br/>")?;
    Ok(())
}

fn range(fh: &mut impl Write,
         color: impl Fn(f64) -> Rgb, n: usize,
         width: u32, comment: &str) -> Result<(), Err> {
    let dt = 1. / (n - 1) as f64;
    let colors: Vec<_> = (0 .. n).map(move |i| color(i as f64 * dt))
        .collect();
    table_of_colors(fh, &colors, width, comment)
}

fn gradient(fh: &mut impl Write, c0: [u8; 3], c1: [u8; 3], n: usize,
            width: u32, comment: &str) -> Result<(), Err> {
    let g = Gradient::new(&Rgb::from(c0), &Rgb::from(c1));
    range(fh, |t| g.rgb(t), n, width, comment)
}

fn palette(fh: &mut impl Write, palette: &Palette, n: usize,
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
    range(&mut fh, |t| Hue().rgb(t), 10, 43, "")?;
    range(&mut fh, |t| Hue().rgb(t), 30, 13, "")?;
    range(&mut fh, |t| Hue().rgb(t), 150, 1, "")?;

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
    for (comment, p) in [("viridis", VIRIDIS.deref()),
                         ("magma", MAGMA.deref()),
                         ("inferno", INFERNO.deref()),
                         ("plasma", PLASMA.deref())] {
        palette(&mut fh, p, 256, 1, comment)?;
    }

    macro_rules! palette { ($($name: ident),*) => {
        $( if let Some(p) = color_brewery::$name.iter().last() {
            let c = format!("{} ({} colors)", stringify!($name), p.len());
            table_of_colors(&mut fh, p.rgb(), 40, &c)?;
            palette(&mut fh, p, 128, 1,
                    &format!("{} (interpolated)", stringify!($name)))?;
        } )*
    }}
    palette!(YLGN, YLGNBU, GNBU, BUGN, PUBUGN, PUBU, BUPU, RDPU, PURD,
             ORRD, YLORRD, YLORBR, PURPLES, BLUES, GREENS, ORANGES,
             REDS, GREYS,
             // Diverging
             PUOR, BRBG, PRGN, PIYG, RDBU, RDGY, RDYLBU, SPECTRAL, RDYLGN);

    writeln!(fh, "<h3>Palettes (others)</h3>")?;
    macro_rules! palette { ($($name: ident),*) => {
        $( if let Some(p) = color_brewery::$name.iter().last() {
            let c = format!("{} ({} colors)", stringify!($name), p.len());
            table_of_colors(&mut fh, p.rgb(), 40, &c)?;
        } )*
    }}
    palette!(SET1, PASTEL1, SET2, PASTEL2, DARK2, SET3, PAIRED, ACCENT);

    writeln!(fh, "</body>\n\
                  </html>")?;
    Ok(())
}
