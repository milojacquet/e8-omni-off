#![allow(unused)]
use crate::e8::Mirror;
use crate::e8::MirrorSet;
use crate::e8::Ring::XX;
use clap::Parser;
use nalgebra::RowSVector;
use nalgebra::SMatrix;
use nalgebra::matrix;
use rand::distr::StandardUniform;
use rand::prelude::*;
use std::ops::Mul;
use std::path::PathBuf;

mod combs;
mod e8;
mod point;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Mirrors
    mirrors: Vec<String>,

    /// Vertex orbits
    #[arg(short, long, value_name = "FILE")]
    vertices: bool,

    /// Write .off
    #[arg(short, long, value_name = "FILE")]
    off: Option<PathBuf>,
}

fn main() {
    let cli = Cli::parse();
    let mut mirror_set = MirrorSet::empty();
    for mirror in cli.mirrors {
        let mirror = mirror.parse().unwrap();
        mirror_set.set(mirror, XX);
    }

    if cli.vertices {
        let vertices = mirror_set.vertex_orbits();
        for vertex in vertices {
            let mut arr: Vec<_> = vertex.rep.iter().copied().collect();
            arr[0] *= vertex.sign;
            println!("{:?}", arr,)
        }
    } else if let Some(file) = cli.off {
        println!("not yet");
    } else {
        println!("what");
    }
}
