#![allow(unused)]
use crate::e8::MirrorSet;
use nalgebra::RowSVector;
use nalgebra::SMatrix;
use nalgebra::matrix;
use rand::distr::StandardUniform;
use rand::prelude::*;
use std::ops::Mul;

mod combs;
mod e8;
mod point;

fn main() {
    print!("{:?}", MirrorSet::e8().vertex())
}
