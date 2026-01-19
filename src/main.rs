#![allow(dead_code)]
use crate::e8::Mirror;
use crate::e8::MirrorSet;
use crate::e8::Ring::XX;
use clap::Parser;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::PathBuf;

mod combs;
mod e8;
mod off;
mod point;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Mirrors
    mirrors: Vec<String>,

    /// Single vertex
    #[arg(short, long)]
    single_vertex: bool,

    /// Vertex orbits
    #[arg(short, long)]
    vertices: bool,

    /// Write .off
    #[arg(short, long, value_name = "FILE")]
    off: Option<PathBuf>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    let mut mirror_set = MirrorSet::empty();
    for mirror in cli.mirrors {
        let mirror = mirror.parse().unwrap();
        mirror_set.set_mirror(mirror, XX);
    }

    if cli.single_vertex {
        let vertex = mirror_set.vertex();
        let mut arr: Vec<_> = vertex.orbit.rep.iter().copied().collect();
        arr[0] *= vertex.orbit.sign;
        println!("{:?}", arr)
    } else if cli.vertices {
        let vertices = mirror_set.vertex_orbits();
        for (vertex, _) in vertices {
            let mut arr: Vec<_> = vertex.orbit.rep.iter().copied().collect();
            arr[0] *= vertex.orbit.sign;
            println!("{:?}", arr)
        }
    } else if let Some(file) = cli.off {
        let mut writer = BufWriter::new(File::create(file)?);
        mirror_set.write_off(&mut writer)?;
        writer.flush()?;
    } else {
        println!("what");
    }
    Ok(())
}
