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
        "use std::sync::LazyLock;\n#[rustfmt::skip]\npub const COMBS:LazyLock<[Vec<Vec<u8>>;9]>=LazyLock::new(|| [{}]);\n",
        combs
            .iter()
            .map(|combs| format!(
                "vec![{}]",
                combs
                    .iter()
                    .map(|combs| format!(
                        "vec![{}]",
                        combs
                            .iter()
                            .map(|comb| comb.to_string())
                            .collect::<Vec<_>>()
                            .join(",")
                    ))
                    .collect::<Vec<_>>()
                    .join(",")
            ))
            .collect::<Vec<_>>()
            .join(",")
    )
    .unwrap();
}
