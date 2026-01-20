use std::array::from_fn;
use std::fs::File;
use std::io::Write;

fn main() {
    let combs: Vec<Vec<Vec<u32>>> = (0..=8)
        .map(|n| {
            (0..=n)
                .map(|k| {
                    (0u32..=(1 << n) - 1)
                        .filter(|i| i.count_ones() == k)
                        .collect()
                })
                .collect()
        })
        .collect();

    let mut file = File::create(concat!(env!("CARGO_MANIFEST_DIR"), "/src/combs.rs")).unwrap();
    write!(
        file,
        "#[rustfmt::skip]
pub const COMBS_LENS:[[u64;9];9]={:?};
#[rustfmt::skip]
pub const COMBS_LISTS:[[[u8;70];9];9]={:?};
#[rustfmt::skip]
pub const COMBS_INDS:[[u64;256];9]={:?};
",
        from_fn::<_, 9, _>(|n| from_fn::<_, 9, _>(|k| combs[n]
            .get(k)
            .map_or(0, |combs| combs.len()))),
        from_fn::<_, 9, _>(|n| from_fn::<_, 9, _>(|k| from_fn::<_, 70, _>(|i| combs[n]
            .get(k)
            .map_or(0, |combs| *combs.get(i).unwrap_or(&0))))),
        from_fn::<_, 9, _>(|n| from_fn::<_, 256, _>(|c| combs[n]
            .get(c.count_ones() as usize)
            .map_or(0, |combs| combs.binary_search(&(c as u32)).unwrap_or(0)))),
    )
    .unwrap();
}
