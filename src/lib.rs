//! Color schemes and gradients.
//!
//! - [`Gradient`]
//! - [`Hue`]
//!
//! [`ColorRange`]
//!
//! Many [`Palette`]s are provided: Matplotlib [`struct@MAGMA`],
//! [`struct@INFERNO`], [`struct@PLASMA`] and [`struct@VIRIDIS`], and
//! all [schemes by Cynthia Brewer](https://colorbrewer2.org/).

use std::f64::consts::PI;
use colorsys::ColorAlpha;

use colorsys::Rgb;

mod palettes;

/// A “continuous” range of colors parametrized by reals in \[0, 1\].
pub trait ColorRange {
    /// Returns the color corresponding to `t` ∈ \[0., 1.\].
    fn rgb(&self, t: f64) -> Rgb;

    /// Return an iterator yielding uniform sampling of `n` points
    /// between `a` and `b` (with the bounds `a` and `b` included in
    /// the list of points) together with colors.  It is not required
    /// that `a <= b`.
    fn range(self, mut a: f64, mut b: f64, n: usize) -> Range<Self>
    where Self: Sized {
        if a == f64::INFINITY { a = f64::MAX; }
        else if a == f64::NEG_INFINITY { a = f64::MIN };
        if b == f64::NEG_INFINITY { b = f64::MIN; }
        else if b == f64::INFINITY { b = f64::MAX };
        // `a` or `b` NaN will give an iterator yielding NaN.
        if n == 0 {
            Range { range: self, a, b, flast: 0., last: 0,
                    i: 1, j: 0 } // Empty iterator
        } else {
            Range { range: self, a, b, flast: (n - 1) as f64,
                    last: n - 1, i: 0, j: n - 1 }
        }
    }
}

/// An iterator yielding `f64` in a given range together with colors.
pub struct Range<R> {
    range: R,
    a: f64, // finite or NaN
    b: f64, // finite if NaN
    flast: f64,
    last: usize,
    i: usize, // first position to be consumed (i ≤ j)
    j: usize, // last position to be consumed
}

impl<R> Range<R> where R: ColorRange {
    /// Return the float and RGB color of the position `k` (assuming
    /// it is in the range `0 ..= self.last`).
    fn rgb(&self, k: usize) -> (f64, Rgb) {
        if k == 0 {
            (self.a, <R as ColorRange>::rgb(&self.range, 0.))
        } else if k == self.last {
            (self.b, <R as ColorRange>::rgb(&self.range, 1.))
        } else {
            let alpha = (self.last - k) as f64;
            let beta = k as f64;
            let t = beta / self.flast;
            let mut x = (alpha * self.a + beta * self.b) / self.flast;
            if x.is_infinite() {
                x = (1. - t) * self.a + t * self.b;
            }
            (x, <R as ColorRange>::rgb(&self.range, t))
        }

    }
}

impl<R> Iterator for Range<R>
where R: ColorRange {
    type Item = (f64, Rgb);

    fn next(&mut self) -> Option<Self::Item> {
        if self.i <= self.j {
            let item = self.rgb(self.i);
            self.i += 1;
            Some(item)
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.last - self.i + 1;
        (len, Some(len))
    }
}

impl<R> ExactSizeIterator for Range<R>
where R: ColorRange {
    fn len(&self) -> usize { self.last - self.i + 1 }
}

impl<R> DoubleEndedIterator for Range<R>
where R: ColorRange {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.i <= self.j {
            let item = self.rgb(self.j);
            if self.j == 0 {
                self.i = 1
            } else {
                self.j -= 1;
            }
            Some(item)
        } else {
            None
        }
    }
}



/// # “Continuous” color ranges
///
/// Color range based on the hue.
///
/// # Example
///
/// ```
/// use color_brewery::{Hue, ColorRange};
/// let rgb = Hue().rgb(0.5);
/// ```
pub struct Hue();

impl ColorRange for Hue {
    /// Return the color corresponding to the hue `h` ∈ \[0., 1.\].
    fn rgb(&self, t: f64) -> Rgb {
        let t = 6. * t;
        let f = 255. * t.fract();
        let ti = t.trunc().rem_euclid(6.);
        if ti == 0.      { Rgb::from((255.,     f,        0.)) }
        else if ti == 1. { Rgb::from((255. - f, 255.,     0.)) }
        else if ti == 2. { Rgb::from((0.,       255.,     f )) }
        else if ti == 3. { Rgb::from((0.,       255. - f, 255.)) }
        else if ti == 4. { Rgb::from((f,        0.,       255.)) }
        else             { Rgb::from((255.,     0.,       255. - f)) }
    }
}

/// The type for colors in the CIE L*C*h*_ab color space with a D50
/// reference white point and an alpha component.  This color space is
/// CIE L*a*b* with polar coordinates.
struct Lch {
    /// The lightness in the range 0. to 100.
    l: f64,
    /// The chroma, in the range 0. to 181.02, but less in practice.
    c: f64,
    /// The hue in degrees in the range 0. to 2π.
    h: f64,
    /// Transparency value ∈ \[0, 1\].
    alpha: f64,
}

const EPS0: f64 = 6. / 29.;
const EPS: f64 = EPS0 * EPS0 * EPS0 ;
const TWO_PI: f64 = 2. * PI;

impl Lch {
    fn from_rgb(c: &Rgb) -> Lch {
        // See https://github.com/dbuenzli/gg/blob/b8704687d669d139bb4ac7a54115afc7e5caaa55/src/gg.ml#L2926
        const C0: f64 = 1. / 3.;
        const C1: f64 = 841. / 108.;
        const C2: f64 = 4. / 29.;
        let r = c.red();
        let g = c.green();
        let b = c.blue();
        let xr = 0.4522795 * r + 0.3993744 * g + 0.1483460 * b;
        let yr = 0.2225105 * r + 0.7168863 * g + 0.0606032 * b;
        let zr = 0.0168820 * r + 0.1176865 * g + 0.8654315 * b;
        let fx = if xr > EPS { xr.powf(C0) } else { C1 * xr + C2 };
        let fy = if yr > EPS { yr.powf(C0) } else { C1 * yr + C2 };
        let fz = if zr > EPS { zr.powf(C0) } else { C1 * zr + C2 };
        let l = 116. * fy - 16.;
        let a = 500. * (fx - fy);
        let b = 200. * (fy - fz);
        let h = { let h = b.atan2(a);
                  if h < 0. { h + TWO_PI } else { h } };
        Lch { l, c: a.hypot(b), h, alpha: c.alpha() }
    }

    fn to_rgb(&self) -> Rgb {
        const C0: f64 = 108. / 841.;
        const C1: f64 = 4. / 29.;
        let a = self.c * self.h.cos();
        let b =  self.c * self.h.sin();
        let fy = (self.l + 16.) / 116.;
        let fx = a / 500. + fy;
        let fz = fy - b / 200.;
        let fx1 = if fx > EPS0 { fx * fx * fx } else { C0 * (fx - C1) };
        let fy1 = if fy > EPS0 { fy * fy * fy } else { C0 * (fy - C1) };
        let fz1 = if fz > EPS0 { fz * fz * fz } else { C0 * (fz - C1) };
        Rgb::new(3.0215932  * fx1 - 1.6168777 * fy1 - 0.4047152 * fz1,
                 -0.9437222 * fx1 + 1.9161365 * fy1 + 0.0275856 * fz1,
                 0.0693906  * fx1 - 0.2290271 * fy1 + 1.1596365 * fz1,
                 Some(self.alpha))
    }
}

/// Gradient between two colors.
///
/// See the [`ColorRange`] trait for methods.
///
/// # Example
///
/// ```
/// use color_brewery::{Rgb, Gradient, ColorRange};
/// let red = Rgb::from((255., 0., 0.));
/// let blue = Rgb::from((0., 0., 255.));
/// let grad = Gradient::new(red, blue);
/// let rgb = grad.rgb(0.5);
/// ```
pub struct Gradient {
    l0: f64, c0: f64,  h0: f64,  a0: f64,
    dl: f64,  dc: f64,  dh: f64,  da: f64,
}

impl Gradient {
    /// Return a gradient from color `c0` to color `c1`.
    pub fn new(c0: &Rgb, c1: &Rgb) -> Gradient {
        let lch0 = Lch::from_rgb(c0);
        let lch1 = Lch::from_rgb(c1);
        let h0 = lch0.h;
        let h1 = lch1.h;
        let dh = {
            if h1 > h0 && h1 - h0 > PI { h1 - (h0 + TWO_PI) }
            else if h1 < h0 && h0 - h1 > PI { h1 + TWO_PI - h0 }
            else { h1 - h0 } };
        Self { l0: lch0.l, c0: lch0.c, h0, a0: lch0.alpha,
               dl: lch1.l - lch0.l,
               dh,
               dc: lch1.c - lch0.c,
               da: lch1.alpha - lch0.alpha }
    }

    /// Returns the color corresponding to `t` ∈ \[0., 1.\] but does
    /// not check the later condition.
    fn rgb_unsafe(&self, t: f64) -> Rgb {
        Lch { l: self.l0 + t * self.dl,
              c: self.c0 + t * self.dc,
              h: self.h0 + t * self.dh,
              alpha: self.a0 + t * self.da
        }.to_rgb()
    }
}

impl ColorRange for Gradient {
    /// Returns the color corresponding to `t` ∈ \[0., 1.\], where
    /// `t == 0.` returns the first color provided in the gradient and
    /// `t == 1.` the second.
    fn rgb(&self, t: f64) -> Rgb { self.rgb_unsafe(t.clamp(0., 1.)) }
}


/// # Color palettes (aka colormaps)
///
/// A Colormap with certain characteristics.
pub use palettes::ty::Palette;

pub use palettes::ty::PaletteType;

impl Palette {
    /// Returns the number of colors in the palette.
    ///
    /// Palettes countains at least 2 colors.
    pub fn len(&self) -> usize { self.rgb.len() }

    /// Says whether the palette is `Seq`uential, `Div`ergent or
    /// `Qual`itative.
    pub fn typ(&self) -> PaletteType { self.typ }

    /// Says whether the palette is colorblind safe (if the palette
    /// specifies it).
    pub fn blind(&self) -> Trivalent { self.blind }

    /// Says whether the palette is suitable for desktop color
    /// printing.
    pub fn print(&self) -> Trivalent { self.print }

    /// Says whether the palette will withstand black and white
    /// photocopying.
    pub fn photocopy(&self) -> Trivalent { self.photocopy }

    /// Says whether the palette is friendly for LCD screens (which
    /// tend to wash-out colors).
    pub fn lcd(&self) -> Trivalent { self.lcd }

    /// Returns the RGB color range of the palette.
    pub fn rgb(& self) -> &Vec<Rgb> { &self.rgb }

    /// Returns a gradient constructed from the palette.
    /// It only makes sense for sequential and some diverging palettes.
    pub fn gradient(&self) -> PaletteGradient {
        PaletteGradient {
            gradients: self.rgb.windows(2)
                .map(|c| Gradient::new(&c[0], &c[1]))
                .collect() }
    }
}

/// A gradient based on a [`Palette`].
pub struct PaletteGradient {
    gradients: Vec<Gradient>,
}

impl ColorRange for PaletteGradient {
    fn rgb(&self, t: f64) -> Rgb {
        let n = self.gradients.len();
        let tn = t.clamp(0., 1.) * n as f64;
        let i = tn.trunc() as usize;
        if i < n { self.gradients[i].rgb_unsafe(tn.fract()) }
        else { self.gradients[n-1].rgb_unsafe(1.) }
    }
}

/// Find palettes matching certain criteria.
pub fn find(len: usize) -> PaletteFind {
    PaletteFind {
        len,
        typ: vec![],
        blind: Trivalent::No, // "no" means "not necessarily want"
        print: Trivalent::No,
        photocopy: Trivalent::No,
        lcd: Trivalent::No,
    }
}

/// Set criteria to find matching palettes.
///
/// Created by [`find`].
#[derive(Clone)]
pub struct PaletteFind {
    len: usize,
    typ: Vec<PaletteType>,
    blind: Trivalent,
    print: Trivalent,
    photocopy: Trivalent,
    lcd: Trivalent,
}

fn satisfy(prop: Trivalent, specified: Trivalent) -> bool {
    use Trivalent::*;
    match specified {
        Yes => matches!(prop, Yes),
        No => true,
        Maybe => matches!(prop, Yes | Maybe),
    }
}

impl PaletteFind {
    /// Find [`Palette`]s with this type.  Use several times to
    /// specify more than one [`PaletteType`].
    pub fn typ(mut self, t: PaletteType) -> Self {
        self.typ.push(t);
        self
    }

    /// Search palettes possibly ([`Maybe`]) or definitely ([`Yes])
    /// suitable for color blind people.
    pub fn blind(mut self, at_least: Trivalent) -> Self {
        self.blind = at_least;
        self
    }

    /// Search palettes possibly ([`Maybe`]) or definitely ([`Yes])
    /// suitable for desktop color printing.
    pub fn print(mut self, at_least: Trivalent) -> Self {
        self.print = at_least;
        self
    }

    /// Search palettes possibly ([`Maybe`]) or definitely ([`Yes])
    /// suitable for black and white photocopying.
    pub fn photocopy(mut self, at_least: Trivalent) -> Self {
        self.photocopy = at_least;
        self
    }

    /// Search palettes possibly ([`Maybe`]) or definitely ([`Yes])
    /// suitable for LCD screens.
    pub fn lcd(mut self, at_least: Trivalent) -> Self {
        self.lcd = at_least;
        self
    }

    /// Return an iterator on all known palettes.
    fn all_palettes() -> impl Iterator<Item = &'static Palette> {
        palettes::ALL_MATPLOTLIB_PALETTES.iter().copied()
            .chain(palettes::ALL_BREWER_PALETTES.iter()
                   .flat_map(|v| v.iter()))
    }

    /// Return the list of palettes of length at least `len` (and
    /// satisfying the criteria set with other methods).
    pub fn palettes(self) -> Vec<&'static Palette> {
        use PaletteType::*;
        let typ = { if self.typ.is_empty() { vec![Seq, Div, Qual] }
                    else { self.typ } };
        Self::all_palettes().filter(|p| {
            p.len() >= self.len
                && typ.contains(&p.typ)
                && satisfy(p.blind, self.blind)
                && satisfy(p.print, self.print)
                && satisfy(p.photocopy, self.photocopy)
                && satisfy(p.lcd, self.lcd)
        }).collect()
    }
}


/// Matplotlib magma color scheme.
///
/// ![magma](https://matplotlib.org/stable/_images/sphx_glr_colormap_reference_001_2_0x.png)
///
pub use palettes::MAGMA;

/// Matplotlib inferno color scheme.
pub use palettes::INFERNO;

/// Matplotlib plasma color scheme.
pub use palettes::PLASMA;

/// Matplotlib viridis color scheme.
pub use palettes::VIRIDIS;


/// Brewer "Light yellow to dark green" sequential scheme
pub use palettes::YLGN;

/// Brewer "Light yellow to green to dark blue" sequential scheme
pub use palettes::YLGNBU;

/// Brewer "Light green to dark blue" sequential scheme
pub use palettes::GNBU;

/// Brewer "Light blue to dark green" sequential scheme
pub use palettes::BUGN;

/// Brewer "Light purple to blue to dark green" sequential scheme
pub use palettes::PUBUGN;

/// Brewer "Light purple to dark blue" sequential scheme
pub use palettes::PUBU;

/// Brewer "Light blue to dark purple" sequential scheme
pub use palettes::BUPU;

/// Brewer "Light red to dark purple" sequential scheme
pub use palettes::RDPU;

/// Brewer "Light purple to dark red" sequential scheme
pub use palettes::PURD;

/// Brewer "Light orange to dark red" sequential scheme
pub use palettes::ORRD;

/// Brewer "Light yellow to orange to dark red" sequential scheme
pub use palettes::YLORRD;

/// Brewer "Light yellow to orange to dark brown" sequential scheme
pub use palettes::YLORBR;


/// Brewer "Light to dark purple" sequential scheme
pub use palettes::PURPLES;

/// Brewer "Light to dark blue" sequential scheme
pub use palettes::BLUES;

/// Brewer "Light to dark green" sequential scheme
pub use palettes::GREENS;

/// Brewer "Light to dark orange" sequential scheme
pub use palettes::ORANGES;

/// Brewer "Light to dark red" sequential scheme
pub use palettes::REDS;

/// Brewer "Light to dark gray" sequential scheme
pub use palettes::GREYS;


/// Brewer "Dark orange to light to dark purple" diverging scheme
pub use palettes::PUOR;

/// Brewer "Dark brown to light to dark blue-green" diverging scheme
pub use palettes::BRBG;

/// Brewer "Dark reddish-purple to light to dark green" diverging scheme
pub use palettes::PRGN;

/// Brewer "Dark magenta to light to dark yellow-green" diverging scheme
pub use palettes::PIYG;

/// Brewer "Dark red to light to dark blue" diverging scheme
pub use palettes::RDBU;

/// Brewer "Dark red to light to dark grey" diverging scheme
pub use palettes::RDGY;

/// Brewer "Dark red to light yelow to dark blue" diverging scheme
pub use palettes::RDYLBU;

/// Brewer "Dark red, orange, light yellow, green, dark blue" diverging scheme
pub use palettes::SPECTRAL;

/// Brewer "Dark red, orange, light yellow, yellow-green, dark green"
/// diverging scheme
pub use palettes::RDYLGN;


/// Brewer qualitative scheme: includes bold, readily named, basic
/// colors (such as red, green, blue).
pub use palettes::SET1;

/// Brewer qualitative scheme: Lighter version of [`struct@SET1`].
pub use palettes::PASTEL1;

/// Brewer qualitative scheme: Includes mostly a mixture colors (such
/// as blue-green, red-orange).
pub use palettes::SET2;

/// Brewer qualitative scheme: Lighter version of [`struct@SET2`].
pub use palettes::PASTEL2;

/// Brewer qualitative scheme: Darker version of [`struct@SET2`].
pub use palettes::DARK2;

/// Brewer qualitative scheme: Medium saturation set with more
/// lightness variation and more classes than [`struct@SET1`] and
/// [`struct@SET2`].
pub use palettes::SET3;

/// Brewer qualitative scheme: Light/dark paris for namable hues.
pub use palettes::PAIRED;

/// Brewer qualitative scheme: Include lightness and saturation
/// extremes to accent small or important areas.
pub use palettes::ACCENT;
use palettes::ty::Trivalent;





// color Blindness
// http://vision.psychol.cam.ac.uk/jdmollon/papers/colourmaps.pdf
// https://www.mapbox.com/blog/colorblind-simulation/
// http://colororacle.org/ — https://github.com/nvkelso/colora-oracle-java



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hue_range() {
        for (i, (x, c)) in Hue().range(0., 1., 11).enumerate() {
            assert!((x - 0.1 * i as f64).abs() <= 1e-15,
                    "{} ≉ {}", x, 0.1 * i as f64);
            assert_eq!(Hue().rgb(x), c);
        }
    }
}
