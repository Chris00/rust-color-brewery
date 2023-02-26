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
use std::marker::PhantomData;
use rgb::{RGBA, RGB8, RGB16, RGBA8, RGBA16};

mod palettes;
use palettes::ty::PaletteData;
pub use palettes::ty::{PaletteType, Trivalent};

/// A “continuous” range of colors parametrized by reals in \[0, 1\].
pub trait ColorRange<Color> {
    /// Returns the color corresponding to `t` ∈ \[0., 1.\].
    fn rgb(&self, t: f64) -> Color;

    /// Return an iterator yielding uniform sampling of `n` points
    /// between `a` and `b` (with the bounds `a` and `b` included in
    /// the list of points) together with colors.  It is not required
    /// that `a <= b`.
    fn range(self, mut a: f64, mut b: f64, n: usize) -> Range<Self, Color>
    where Self: Sized {
        if a == f64::INFINITY { a = f64::MAX; }
        else if a == f64::NEG_INFINITY { a = f64::MIN };
        if b == f64::NEG_INFINITY { b = f64::MIN; }
        else if b == f64::INFINITY { b = f64::MAX };
        // `a` or `b` NaN will give an iterator yielding NaN.
        if n == 0 {
            Range { range: self,  color: PhantomData,
                    a, b, flast: 0., last: 0,
                    i: 1, j: 0 } // Empty iterator
        } else {
            Range { range: self,  color: PhantomData,
                    a, b, flast: (n - 1) as f64,
                    last: n - 1, i: 0, j: n - 1 }
        }
    }
}

/// An iterator yielding `f64` in a given range together with colors.
pub struct Range<R, Color> {
    range: R,
    color: PhantomData<Color>,
    a: f64, // finite or NaN
    b: f64, // finite or NaN
    flast: f64, // `last` as a floating-point number
    last: usize,
    i: usize, // first position to be consumed (i ≤ j)
    j: usize, // last position to be consumed
}

impl<R, Color> Range<R, Color> where R: ColorRange<Color> {
    /// Return the float and RGB color of the position `k` (assuming
    /// it is in the range `0 ..= self.last`).
    fn rgb(&self, k: usize) -> (f64, Color) {
        if k == 0 {
            (self.a, R::rgb(&self.range, 0.))
        } else if k == self.last {
            (self.b, R::rgb(&self.range, 1.))
        } else {
            let alpha = (self.last - k) as f64;
            let beta = k as f64;
            let t = beta / self.flast;
            let mut x = (alpha * self.a + beta * self.b) / self.flast;
            if x.is_infinite() {
                x = (1. - t) * self.a + t * self.b;
            }
            (x, R::rgb(&self.range, t))
        }

    }
}

impl<R, Color> Iterator for Range<R, Color>
where R: ColorRange<Color> {
    type Item = (f64, Color);

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

impl<R, Color> ExactSizeIterator for Range<R, Color>
where R: ColorRange<Color> {
    fn len(&self) -> usize { self.last - self.i + 1 }
}

impl<R, Color> DoubleEndedIterator for Range<R, Color>
where R: ColorRange<Color> {
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

/// Specifies the methods a RGB color encoding must provide.
pub trait RGBColor: Sized {
    /// Return the red, green, blue and alpha components of the color
    /// (in \[0, 255\]).
    fn to_rgba(&self) -> RGBA<f64>;

    /// Create a color from its RGBA components (in \[0, 255\]).
    fn from_rgba(rgba: RGBA<f64>) -> Self;

    /// Return the color corresponding to the hue `h` ∈ \[0., 1.\].
    ///
    /// # Example
    ///
    /// ```
    /// use rgb::RGB8;
    /// use color_brewery::{RGBColor, ColorRange};
    /// let rgb = RGB8::HUE.rgb(0.5);
    /// ```
    const HUE: Hue<Self> = Hue { color: PhantomData };

    /// Return a gradient from color `c0` to color `c1`.
    ///
    /// # Example
    ///
    /// ```
    /// use rgb::RGB8;
    /// use color_brewery::{RGBColor, ColorRange};
    /// let red = RGB8::new(255,0, 0);
    /// let blue = RGB8::new(0, 0, 255);
    /// let grad = red.gradient(&blue);
    /// let rgb = grad.rgb(0.5);
    /// ```
    fn gradient(&self, c1: &Self) -> Gradient<Self> {
        let lch0 = Lch::from_rgb(Self::to_rgba(self));
        let lch1 = Lch::from_rgb(Self::to_rgba(c1));
        let h0 = lch0.h;
        let h1 = lch1.h;
        let dh = {
            if h1 > h0 && h1 - h0 > PI { h1 - (h0 + TWO_PI) }
            else if h1 < h0 && h0 - h1 > PI { h1 + TWO_PI - h0 }
            else { h1 - h0 } };
        Gradient { c0: lch0,
                   dc: Lch { l: lch1.l - lch0.l, c: lch1.c - lch0.c,
                             h: dh, a: lch1.a - lch0.a },
                   color: PhantomData }
    }

    /// Find palettes matching certain criteria.
    fn palettes(len: usize) -> PaletteFind<Self> {
        PaletteFind {
            len,
            typ: vec![],
            blind: Trivalent::No, // "no" means "not necessarily want"
            print: Trivalent::No,
            photocopy: Trivalent::No,
            lcd: Trivalent::No,
            color: PhantomData,
        }
    }


    /// Matplotlib magma color scheme.
    ///
    /// ![magma](https://matplotlib.org/stable/_images/sphx_glr_colormap_reference_001_2_0x.png)
    #[inline]
    fn magma() -> Palette<Self> { Palette::new(&palettes::MAGMA) }

    /// Matplotlib inferno color scheme.
    #[inline]
    fn inferno() -> Palette<Self> { Palette::new(&palettes::INFERNO) }

    /// Matplotlib plasma color scheme.
    #[inline]
    fn plasma() -> Palette<Self> { Palette::new(&palettes::PLASMA) }

    /// Matplotlib viridis color scheme.
    #[inline]
    fn viridis() -> Palette<Self> { Palette::new(&palettes::VIRIDIS) }

    /// Brewer "Light yellow to dark green" sequential scheme.
    #[inline]
    fn ylgn() -> PaletteIter<Self> { PaletteIter::new(&palettes::YLGN) }

    /// Brewer "Light yellow to green to dark blue" sequential scheme.
    #[inline]
    fn ylgnbu() -> PaletteIter<Self> { PaletteIter::new(&palettes::YLGNBU) }

    /// Brewer "Light green to dark blue" sequential scheme.
    #[inline]
    fn gnbu() -> PaletteIter<Self> { PaletteIter::new(&palettes::GNBU) }

    /// Brewer "Light blue to dark green" sequential scheme.
    #[inline]
    fn bugn() -> PaletteIter<Self> { PaletteIter::new(&palettes::BUGN) }

    /// Brewer "Light purple to blue to dark green" sequential scheme.
    #[inline]
    fn pubugn() -> PaletteIter<Self> { PaletteIter::new(&palettes::PUBUGN) }

    /// Brewer "Light purple to dark blue" sequential scheme.
    #[inline]
    fn pubu() -> PaletteIter<Self> { PaletteIter::new(&palettes::PUBU) }

    /// Brewer "Light blue to dark purple" sequential scheme.
    #[inline]
    fn bupu() -> PaletteIter<Self> { PaletteIter::new(&palettes::BUPU) }

    /// Brewer "Light red to dark purple" sequential scheme.
    #[inline]
    fn rdpu() -> PaletteIter<Self> { PaletteIter::new(&palettes::RDPU) }

    /// Brewer "Light purple to dark red" sequential scheme.
    #[inline]
    fn purd() -> PaletteIter<Self> { PaletteIter::new(&palettes::PURD) }

    /// Brewer "Light orange to dark red" sequential scheme.
    #[inline]
    fn orrd() -> PaletteIter<Self> { PaletteIter::new(&palettes::ORRD) }

    /// Brewer "Light yellow to orange to dark red" sequential scheme.
    #[inline]
    fn ylorrd() -> PaletteIter<Self> { PaletteIter::new(&palettes::YLORRD) }

    /// Brewer "Light yellow to orange to dark brown" sequential scheme.
    #[inline]
    fn ylorbr() -> PaletteIter<Self> { PaletteIter::new(&palettes::YLORBR) }


    /// Brewer "Light to dark purple" sequential scheme.
    #[inline]
    fn purples() -> PaletteIter<Self> { PaletteIter::new(&palettes::PURPLES) }

    /// Brewer "Light to dark blue" sequential scheme.
    #[inline]
    fn blues() -> PaletteIter<Self> { PaletteIter::new(&palettes::BLUES) }

    /// Brewer "Light to dark green" sequential scheme.
    #[inline]
    fn greens() -> PaletteIter<Self> { PaletteIter::new(&palettes::GREENS) }

    /// Brewer "Light to dark orange" sequential scheme.
    #[inline]
    fn oranges() -> PaletteIter<Self> { PaletteIter::new(&palettes::ORANGES) }

    /// Brewer "Light to dark red" sequential scheme.
    #[inline]
    fn reds() -> PaletteIter<Self> { PaletteIter::new(&palettes::REDS) }

    /// Brewer "Light to dark gray" sequential scheme.
    #[inline]
    fn greys() -> PaletteIter<Self> { PaletteIter::new(&palettes::GREYS) }


    /// Brewer "Dark orange to light to dark purple" diverging scheme
    #[inline]
    fn puor() -> PaletteIter<Self> { PaletteIter::new(&palettes::PUOR) }

    /// Brewer "Dark brown to light to dark blue-green" diverging scheme
    #[inline]
    fn brbg() -> PaletteIter<Self> { PaletteIter::new(&palettes::BRBG) }

    /// Brewer "Dark reddish-purple to light to dark green" diverging scheme
    #[inline]
    fn prgn() -> PaletteIter<Self> { PaletteIter::new(&palettes::PRGN) }

    /// Brewer "Dark magenta to light to dark yellow-green" diverging scheme
    #[inline]
    fn piyg() -> PaletteIter<Self> { PaletteIter::new(&palettes::PIYG) }

    /// Brewer "Dark red to light to dark blue" diverging scheme.
    #[inline]
    fn rdbu() -> PaletteIter<Self> { PaletteIter::new(&palettes::RDBU) }

    /// Brewer "Dark red to light to dark grey" diverging scheme.
    #[inline]
    fn rdgy() -> PaletteIter<Self> { PaletteIter::new(&palettes::RDGY) }

    /// Brewer "Dark red to light yelow to dark blue" diverging scheme.
    #[inline]
    fn rdylbu() -> PaletteIter<Self> { PaletteIter::new(&palettes::RDYLBU) }

    /// Brewer "Dark red, orange, light yellow, green, dark blue"
    /// diverging scheme.
    #[inline]
    fn spectral() -> PaletteIter<Self> { PaletteIter::new(&palettes::SPECTRAL) }

    /// Brewer "Dark red, orange, light yellow, yellow-green, dark green"
    /// diverging scheme.
    #[inline]
    fn rdylgn() -> PaletteIter<Self> { PaletteIter::new(&palettes::RDYLGN) }


    /// Brewer qualitative scheme: includes bold, readily named, basic
    /// colors (such as red, green, blue).
    #[inline]
    fn set1() -> PaletteIter<Self> { PaletteIter::new(&palettes::SET1) }

    /// Brewer qualitative scheme: Lighter version of [`struct@SET1`].
    #[inline]
    fn pastel1() -> PaletteIter<Self> { PaletteIter::new(&palettes::PASTEL1) }

    /// Brewer qualitative scheme: Includes mostly a mixture colors
    /// (such as blue-green, red-orange).
    #[inline]
    fn set2() -> PaletteIter<Self> { PaletteIter::new(&palettes::SET2) }

    /// Brewer qualitative scheme: Lighter version of [`struct@SET2`].
    #[inline]
    fn pastel2() -> PaletteIter<Self> { PaletteIter::new(&palettes::PASTEL2) }

    /// Brewer qualitative scheme: Darker version of [`struct@SET2`].
    #[inline]
    fn dark2() -> PaletteIter<Self> { PaletteIter::new(&palettes::DARK2) }

    /// Brewer qualitative scheme: Medium saturation set with more
    /// lightness variation and more classes than [`struct@SET1`] and
    /// [`struct@SET2`].
    #[inline]
    fn set3() -> PaletteIter<Self> { PaletteIter::new(&palettes::SET3) }

    /// Brewer qualitative scheme: Light/dark paris for namable hues.
    #[inline]
    fn paired() -> PaletteIter<Self> { PaletteIter::new(&palettes::PAIRED) }

    /// Brewer qualitative scheme: Include lightness and saturation
    /// extremes to accent small or important areas.
    #[inline]
    fn accent() -> PaletteIter<Self> { PaletteIter::new(&palettes::ACCENT) }

    /// Convert the color to grayscale.
    fn to_gray(&self) -> Self {
        let RGBA{ r, g, b, a } = Self::to_rgba(self);
        let x = 0.299 * r + 0.587 * g + 0.114 * b;
        Self::from_rgba(RGBA{ r: x, g: x, b: x, a })
    }
}

impl RGBColor for RGBA<f64> {
    #[inline]
    fn to_rgba(&self) -> RGBA<f64> { *self }

    #[inline]
    fn from_rgba(c: RGBA<f64>) -> Self { c }
}

impl RGBColor for RGB8 {
    #[inline]
    fn to_rgba(&self) -> RGBA<f64> {
        RGBA{ r: self.r as f64, g: self.g as f64, b: self.b as f64, a: 255. }
    }

    #[inline]
    fn from_rgba(c: RGBA<f64>) -> Self {
        RGB8 { r: c.r as u8,  g: c.g as u8,  b: c.b as u8 }
    }
}

impl RGBColor for RGB16 {
    #[inline]
    fn to_rgba(&self) -> RGBA<f64> {
        RGBA{ r: self.r as f64, g: self.g as f64, b: self.b as f64, a: 255. }
    }

    #[inline]
    fn from_rgba(c: RGBA<f64>) -> Self {
        RGB16 { r: c.r as u16,  g: c.g as u16,  b: c.b as u16 }
    }
}

impl RGBColor for RGBA8 {
    #[inline]
    fn to_rgba(&self) -> RGBA<f64> {
        RGBA{ r: self.r as f64, g: self.g as f64, b: self.b as f64,
              a: self.a as f64 }
    }

    #[inline]
    fn from_rgba(c: RGBA<f64>) -> Self {
        RGBA8 { r: c.r as u8,  g: c.g as u8,  b: c.b as u8, a: c.a as u8 }
    }
}

impl RGBColor for RGBA16 {
    #[inline]
    fn to_rgba(&self) -> RGBA<f64> {
        RGBA{ r: self.r as f64, g: self.g as f64, b: self.b as f64,
              a: self.a as f64 }
    }

    #[inline]
    fn from_rgba(c: RGBA<f64>) -> Self {
        RGBA16 { r: c.r as u16,  g: c.g as u16,  b: c.b as u16, a: c.a as u16 }
    }
}

/// The type for colors in the CIE L*C*h*_ab color space with a D50
/// reference white point and an alpha component.  This color space is
/// CIE L*a*b* with polar coordinates.
#[derive(Clone, Copy)]
struct Lch {
    /// The lightness in the range 0. to 100.
    l: f64,
    /// The chroma, in the range 0. to 181.02, but less in practice.
    c: f64,
    /// The hue in degrees in the range 0. to 2π.
    h: f64,
    /// Alpha component
    a: f64,
}

const EPS0: f64 = 6. / 29.;
const EPS: f64 = EPS0 * EPS0 * EPS0 ;
const TWO_PI: f64 = 2. * PI;

impl Lch {
    fn from_rgb(c: RGBA<f64>) -> Lch {
        // See https://github.com/dbuenzli/gg/blob/b8704687d669d139bb4ac7a54115afc7e5caaa55/src/gg.ml#L2926
        const C0: f64 = 1. / 3.;
        const C1: f64 = 841. / 108.;
        const C2: f64 = 4. / 29.;
        let xr = 0.4522795 * c.r + 0.3993744 * c.g + 0.1483460 * c.b;
        let yr = 0.2225105 * c.r + 0.7168863 * c.g + 0.0606032 * c.b;
        let zr = 0.0168820 * c.r + 0.1176865 * c.g + 0.8654315 * c.b;
        let fx = if xr > EPS { xr.powf(C0) } else { C1 * xr + C2 };
        let fy = if yr > EPS { yr.powf(C0) } else { C1 * yr + C2 };
        let fz = if zr > EPS { zr.powf(C0) } else { C1 * zr + C2 };
        let l = 116. * fy - 16.;
        let a = 500. * (fx - fy);
        let b = 200. * (fy - fz);
        let h = { let h = b.atan2(a);
                  if h < 0. { h + TWO_PI } else { h } };
        Lch { l, c: a.hypot(b), h, a: c.a }
    }

    fn to_rgb(&self) -> RGBA<f64> {
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
        let r = 3.0215932  * fx1 - 1.6168777 * fy1 - 0.4047152 * fz1;
        let g = -0.9437222 * fx1 + 1.9161365 * fy1 + 0.0275856 * fz1;
        let b = 0.0693906  * fx1 - 0.2290271 * fy1 + 1.1596365 * fz1;
        RGBA { r, g, b, a: self.a }
    }
}

/// Hue
///
pub struct Hue<Color> { color: PhantomData<Color> }

impl<Color: RGBColor> ColorRange<Color> for Hue<Color> {
    fn rgb(&self, t: f64) -> Color {
        let t = 6. * t;
        let f = 255. * t.fract();
        let ti = t.trunc().rem_euclid(6.);
        let rgba = {
            if ti == 0.      { RGBA{ r: 255., g: f,     b: 0.,      a: 255.} }
            else if ti == 1. { RGBA{ r: 255. - f, g: 255., b: 0.,   a: 255.} }
            else if ti == 2. { RGBA{ r: 0.,   g: 255.,  b: f,       a: 255.} }
            else if ti == 3. { RGBA{ r: 0.,   g: 255. - f, b: 255., a: 255.} }
            else if ti == 4. { RGBA{ r: f,    g: 0.,    b: 255.,    a: 255.} }
            else             { RGBA{ r: 255., g: 0.,    b: 255. - f, a: 255.} }
        };
        Color::from_rgba(rgba)
    }
}


/// Gradient between two colors.
///
/// Created by [`RGBcolor::gradient`].  See the [`ColorRange`] trait
/// for methods.
pub struct Gradient<Color> {
    c0: Lch, // first color
    dc: Lch, // last - fist color
    color: PhantomData<Color>,
}

impl<Color> Gradient<Color>
where Color: RGBColor {
    /// Returns the color corresponding to `t` ∈ \[0., 1.\] but does
    /// not check the later condition.
    #[inline]
    fn rgb_unsafe(&self, t: f64) -> Color {
        let lhc = Lch { l: self.c0.l + t * self.dc.l,
              c: self.c0.c + t * self.dc.c,
              h: self.c0.h + t * self.dc.h,
              a: self.c0.a + t * self.dc.a };
        Color::from_rgba(lhc.to_rgb())
    }
}

impl<Color> ColorRange<Color> for Gradient<Color>
where Color: RGBColor {
    /// Returns the color corresponding to `t` ∈ \[0., 1.\], where
    /// `t == 0.` returns the first color provided in the gradient and
    /// `t == 1.` the second.
    fn rgb(&self, t: f64) -> Color { self.rgb_unsafe(t.clamp(0., 1.)) }
}


#[derive(Clone, Copy)]
pub struct Palette<Color> {
    palette: &'static PaletteData,
    color: PhantomData<Color>,
}

impl<Color: RGBColor> Palette<Color> {
    fn new(palette: &'static PaletteData) -> Self {
        Self { palette, color: PhantomData }
    }
}

/// # Color palettes (aka colormaps)
///
/// A Colormap with certain characteristics.
impl<Color> Palette<Color>
where Color: RGBColor {
    /// Returns the number of colors in the palette.
    ///
    /// Palettes countains at least 2 colors.
    pub fn len(&self) -> usize { self.palette.rgb.len() }

    /// Says whether the palette is `Seq`uential, `Div`ergent or
    /// `Qual`itative.
    pub fn typ(&self) -> PaletteType { self.palette.typ }

    /// Says whether the palette is colorblind safe (if the palette
    /// specifies it).
    pub fn blind(&self) -> Trivalent { self.palette.blind }

    /// Says whether the palette is suitable for desktop color
    /// printing.
    pub fn print(&self) -> Trivalent { self.palette.print }

    /// Says whether the palette will withstand black and white
    /// photocopying.
    pub fn photocopy(&self) -> Trivalent { self.palette.photocopy }

    /// Says whether the palette is friendly for LCD screens (which
    /// tend to wash-out colors).
    pub fn lcd(&self) -> Trivalent { self.palette.lcd }

    /// Returns the RGB colors of the palette.
    pub fn colors(&self) -> Vec<Color> {
        self.palette.rgb.iter().map(|&c| Color::from_rgba(c)).collect()
    }

    /// Returns a gradient constructed from the palette.
    /// It only makes sense for sequential and some diverging palettes.
    pub fn gradient(&self) -> PaletteGradient<Color> {
        PaletteGradient {
            gradients: self.palette.rgb.windows(2)
                .map(|c| { let c0 = Color::from_rgba(c[0]);
                           let c1 = Color::from_rgba(c[1]);
                           c0.gradient(&c1) })
                .collect() }
    }
}

/// A gradient based on a [`Palette`].
pub struct PaletteGradient<Color> {
    gradients: Vec<Gradient<Color>>,
}

impl<Color> ColorRange<Color> for PaletteGradient<Color>
where Color: RGBColor {
    fn rgb(&self, t: f64) -> Color {
        let n = self.gradients.len();
        let tn = t.clamp(0., 1.) * n as f64;
        let i = tn.trunc() as usize;
        if i < n { self.gradients[i].rgb_unsafe(tn.fract()) }
        else { self.gradients[n-1].rgb_unsafe(1.) }
    }
}

/// An exact size iterator over [`Palette`]s.
#[derive(Clone, Copy)]
pub struct PaletteIter<Color> {
    palettes: &'static Vec<PaletteData>,
    i: usize, // first position to be consumed (i < j)
    j: usize, // position after the last one to be consumed
    color: PhantomData<Color>,
}

impl<Color: RGBColor> PaletteIter<Color> {
    fn new(palettes: &'static Vec<PaletteData>) -> Self {
        Self { palettes, i: 0, j: palettes.len(),
               color: PhantomData }
    }
}

impl<Color> Iterator for PaletteIter<Color>
where Color: RGBColor {
    type Item = Palette<Color>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.i >= self.j { return None }
        let x = Palette::new(&self.palettes[self.i]);
        self.i += 1;
        Some(x)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.j - self.i;
        (len, Some(len))
    }
}

impl<Color: RGBColor> ExactSizeIterator for PaletteIter<Color> {}

impl<Color: RGBColor> DoubleEndedIterator for PaletteIter<Color>{
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.i >= self.j { return None }
        self.j -= 1;
        Some(Palette::new(&self.palettes[self.j]))
    }
}


/// Set criteria to find matching palettes.
///
/// Created by [`find`].
#[derive(Clone)]
pub struct PaletteFind<Color> {
    len: usize,
    typ: Vec<PaletteType>,
    blind: Trivalent,
    print: Trivalent,
    photocopy: Trivalent,
    lcd: Trivalent,
    color: PhantomData<Color>,
}

fn satisfy(prop: Trivalent, specified: Trivalent) -> bool {
    use Trivalent::*;
    match specified {
        Yes => matches!(prop, Yes),
        No => true,
        Maybe => matches!(prop, Yes | Maybe),
    }
}

impl<Color> PaletteFind<Color>
where Color: RGBColor {
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
    fn all() -> impl Iterator<Item = Palette<Color>> {
        palettes::ALL_MATPLOTLIB_PALETTES.iter().copied()
            .chain(palettes::ALL_BREWER_PALETTES.iter()
                   .flat_map(|v| v.iter()))
            .map(|c| Palette::new(c))
    }

    /// Return the list of palettes of length at least `len` (and
    /// satisfying the criteria set with other methods).
    pub fn find(self) -> impl Iterator<Item = Palette<Color>> {
        use PaletteType::*;
        let typ = { if self.typ.is_empty() { vec![Seq, Div, Qual] }
                    else { self.typ } };
        Self::all().filter(move |p| {
            p.len() >= self.len
                && typ.contains(&p.palette.typ)
                && satisfy(p.palette.blind, self.blind)
                && satisfy(p.palette.print, self.print)
                && satisfy(p.palette.photocopy, self.photocopy)
                && satisfy(p.palette.lcd, self.lcd)
        })
    }
}



// color Blindness
// http://vision.psychol.cam.ac.uk/jdmollon/papers/colourmaps.pdf
// https://www.mapbox.com/blog/colorblind-simulation/
// http://colororacle.org/ — https://github.com/nvkelso/colora-oracle-java



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hue_range() {
        for (i, (x, c)) in RGB8::HUE.range(0., 1., 11).enumerate() {
            assert!((x - 0.1 * i as f64).abs() <= 1e-15,
                    "{} ≉ {}", x, 0.1 * i as f64);
            assert_eq!(RGB8::HUE.rgb(x), c);
        }
    }
}
