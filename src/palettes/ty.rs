pub(crate) use rgb::RGBA;

/// A Colormap with certain characteristics.
pub(crate) struct PaletteData {
    pub(crate) rgb: Vec<RGBA<f64>>, // Invariant: length ≥ 2
    pub(crate) typ: PaletteType,
    pub(crate) blind: Trivalent,
    pub(crate) print: Trivalent,
    pub(crate) photocopy: Trivalent,
    pub(crate) lcd: Trivalent,
}

/// Type of Palette.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PaletteType {
    /// Sequential color scheme, suited to ordered data that progress
    /// from low to high. Lightness steps dominate the look of these
    /// schemes, with light colors for low data values to dark colors
    /// for high data values.
    Seq,
    /// Divergent color scheme.  They put equal emphasis on mid-range
    /// critical values and extremes at both ends of the data
    /// range.  The critical class or break in the middle of the legend
    /// is emphasized with light colors and low and high extremes are
    /// emphasized with dark colors that have contrasting hues.
    Div,
    /// Qualitative color scheme.  They do not imply magnitude
    /// differences between legend classes, and hues are used to
    /// create the primary visual differences between classes.
    /// Qualitative schemes are best suited to representing nominal or
    /// categorical data.
    Qual
}

/// Trivalent logic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Trivalent {
    Yes,
    Maybe,
    No,
}
