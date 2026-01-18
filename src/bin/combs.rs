use std::fs::File;
use std::io::Write;

fn main() {
    let combs: Vec<Vec<Vec<u32>>> = (0..=8)
        .map(|n| {
            (0..=n)
                .map(|k| (0u32..=(1 << n) - 1).map(|i| i.count_ones()).collect())
                .collect()
        })
        .collect();

    let combs_str = todo!();

    let mut file = File::create("../combs.rs").unwrap();
    write!(file, "const ").unwrap();
}
