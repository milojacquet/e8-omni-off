use crate::e8::Ring::XX;
use crate::e8::Ring::oo;
use crate::point::Point;
use crate::point::Vec8;
use bitflags::bitflags;
use fxhash::FxHashMap;
use fxhash::FxHasher;
use nalgebra::RowSVector;
use nalgebra::SMatrix;
use nalgebra::matrix;
use rand::distr::StandardUniform;
use rand::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::ops::Mul;
use std::str::FromStr;

const E8_SIZE: u64 = 696729600;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Mirror {
    A0,
    A1,
    A2,
    A3,
    B0,
    B1,
    C,
    M,
}

impl Mirror {
    const ALL: [Self; 8] = [
        Self::A0,
        Self::A1,
        Self::A2,
        Self::A3,
        Self::B0,
        Self::B1,
        Self::C,
        Self::M,
    ];

    /// reflection pole with norm 2√2
    pub fn pole(&self) -> Vec8 {
        match self {
            Mirror::C => matrix![2, -2, 0, 0, 0, 0, 0, 0],
            Mirror::M => matrix![0, 2, -2, 0, 0, 0, 0, 0],
            Mirror::A3 => matrix![0, 0, 2, -2, 0, 0, 0, 0],
            Mirror::A2 => matrix![0, 0, 0, 2, -2, 0, 0, 0],
            Mirror::A1 => matrix![0, 0, 0, 0, 2, -2, 0, 0],
            Mirror::A0 => matrix![0, 0, 0, 0, 0, 2, -2, 0],
            Mirror::B1 => matrix![-2, -2, 0, 0, 0, 0, 0, 0],
            Mirror::B0 => matrix![1, 1, 1, 1, 1, 1, 1, 1],
        }
    }

    pub fn mat(&self) -> E8 {
        E8(SMatrix::identity() * 4 - self.pole().transpose() * self.pole())
    }
}

impl FromStr for Mirror {
    type Err = ();
    fn from_str(st: &str) -> Result<Self, ()> {
        match st {
            "A0" => Ok(Mirror::A0),
            "A1" => Ok(Mirror::A1),
            "A2" => Ok(Mirror::A2),
            "A3" => Ok(Mirror::A3),
            "B0" => Ok(Mirror::B0),
            "B1" => Ok(Mirror::B1),
            "C" => Ok(Mirror::C),
            "M" => Ok(Mirror::M),
            _ => Err(()),
        }
    }
}

/// Matrix in E8 group times 4
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct E8(SMatrix<i16, 8, 8>);

impl E8 {
    pub fn identity() -> Self {
        E8(SMatrix::identity() * 2)
    }

    pub fn inv(self) -> Self {
        E8(self.0.transpose())
    }
}

impl Mul<E8> for E8 {
    type Output = E8;
    fn mul(self, other: E8) -> E8 {
        E8(self.0 * other.0 / 4)
    }
}

impl Mul<E8> for Point {
    type Output = Point;
    fn mul(self, other: E8) -> Point {
        Point::new(self.vec() * other.0 / 4)
    }
}

impl Distribution<E8> for StandardUniform {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> E8 {
        MirrorSet::all().sample(rng)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Ring {
    #[expect(non_camel_case_types)]
    oo,
    XX,
}

impl Ring {
    pub fn size(self) -> u8 {
        match self {
            oo => 0,
            XX => 1,
        }
    }

    pub fn toggle(self) -> Self {
        match self {
            oo => XX,
            XX => oo,
        }
    }
}

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
    pub struct MirrorSet: u8 {
        const A0 = 1 << 0;
        const A1 = 1 << 1;
        const A2 = 1 << 2;
        const A3 = 1 << 3;
        const B0 = 1 << 4;
        const B1 = 1 << 5;
        const C = 1 << 6;
        const M = 1 << 7;
    }
}

impl MirrorSet {
    pub fn a0(self) -> Ring {
        if self.contains(Self::A0) { XX } else { oo }
    }
    pub fn a1(self) -> Ring {
        if self.contains(Self::A1) { XX } else { oo }
    }
    pub fn a2(self) -> Ring {
        if self.contains(Self::A2) { XX } else { oo }
    }
    pub fn a3(self) -> Ring {
        if self.contains(Self::A3) { XX } else { oo }
    }
    pub fn b0(self) -> Ring {
        if self.contains(Self::B0) { XX } else { oo }
    }
    pub fn b1(self) -> Ring {
        if self.contains(Self::B1) { XX } else { oo }
    }
    pub fn c(self) -> Ring {
        if self.contains(Self::C) { XX } else { oo }
    }
    pub fn m(self) -> Ring {
        if self.contains(Self::M) { XX } else { oo }
    }

    pub fn mirror(mirror: Mirror) -> Self {
        match mirror {
            Mirror::A0 => Self::A0,
            Mirror::A1 => Self::A1,
            Mirror::A2 => Self::A2,
            Mirror::A3 => Self::A3,
            Mirror::B0 => Self::B0,
            Mirror::B1 => Self::B1,
            Mirror::C => Self::C,
            Mirror::M => Self::M,
        }
    }

    pub fn has_mirror(self, mirror: Mirror) -> bool {
        match mirror {
            Mirror::A0 => self.a0() == XX,
            Mirror::A1 => self.a1() == XX,
            Mirror::A2 => self.a2() == XX,
            Mirror::A3 => self.a3() == XX,
            Mirror::B0 => self.b0() == XX,
            Mirror::B1 => self.b1() == XX,
            Mirror::C => self.c() == XX,
            Mirror::M => self.m() == XX,
        }
    }

    pub fn set_mirror(&mut self, mirror: Mirror, val: Ring) {
        match mirror {
            Mirror::A0 => self.set(Self::A0, val == XX),
            Mirror::A1 => self.set(Self::A1, val == XX),
            Mirror::A2 => self.set(Self::A2, val == XX),
            Mirror::A3 => self.set(Self::A3, val == XX),
            Mirror::B0 => self.set(Self::B0, val == XX),
            Mirror::B1 => self.set(Self::B1, val == XX),
            Mirror::C => self.set(Self::C, val == XX),
            Mirror::M => self.set(Self::M, val == XX),
        }
    }

    pub fn order(self) -> u64 {
        let mut order = 1;

        let (o, len_a) = match (self.a0(), self.a1(), self.a2(), self.a3(), self.m()) {
            (oo, oo, oo, oo, _) => (1, 0),
            (XX, oo, oo, oo, _) => (2, 0),
            (oo, XX, oo, oo, _) => (2, 0),
            (XX, XX, oo, oo, _) => (6, 0),
            (oo, oo, XX, oo, _) => (2, 0),
            (XX, oo, XX, oo, _) => (2 * 2, 0),
            (oo, XX, XX, oo, _) => (6, 0),
            (XX, XX, XX, oo, _) => (24, 0),
            (oo, oo, oo, XX, oo) => (2, 0),
            (XX, oo, oo, XX, oo) => (2 * 2, 0),
            (oo, XX, oo, XX, oo) => (2 * 2, 0),
            (XX, XX, oo, XX, oo) => (6 * 2, 0),
            (oo, oo, XX, XX, oo) => (6, 0),
            (XX, oo, XX, XX, oo) => (2 * 6, 0),
            (oo, XX, XX, XX, oo) => (24, 0),
            (XX, XX, XX, XX, oo) => (120, 0),
            (oo, oo, oo, XX, XX) => (1, 1),
            (XX, oo, oo, XX, XX) => (2, 1),
            (oo, XX, oo, XX, XX) => (2, 1),
            (XX, XX, oo, XX, XX) => (6, 1),
            (oo, oo, XX, XX, XX) => (1, 2),
            (XX, oo, XX, XX, XX) => (2, 2),
            (oo, XX, XX, XX, XX) => (1, 3),
            (XX, XX, XX, XX, XX) => (1, 4),
        };
        order *= o;

        let (o, len_b) = match (self.b0(), self.b1(), self.m()) {
            (oo, oo, oo) => (1, 0),
            (XX, oo, oo) => (2, 0),
            (oo, XX, oo) => (2, 1),
            (XX, XX, oo) => (6, 2),
            (oo, oo, XX) => (1, 0),
            (XX, oo, XX) => (2, 0),
            (oo, XX, XX) => (1, 1),
            (XX, XX, XX) => (1, 2),
        };
        order *= o;

        let (o, len_c) = match (self.c(), self.m()) {
            (oo, oo) => (1, 0),
            (XX, oo) => (2, 0),
            (oo, XX) => (1, 0),
            (XX, XX) => (1, 1),
        };
        order *= o;

        let mut lens = [len_a, len_b, len_c];
        lens.sort();
        if self.m() == XX {
            order *= match lens {
                [0, 0, 0] => 2,
                [0, 0, 1] => 6,
                [0, 1, 1] => 24,
                [1, 1, 1] => 192,
                [0, 0, 2] => 24,
                [0, 1, 2] => 120,
                [1, 1, 2] => 1920,
                [0, 2, 2] => 720,
                [1, 2, 2] => 51840,
                [0, 0, 3] => 120,
                [0, 1, 3] => 720,
                [1, 1, 3] => 32 * 720,
                [0, 2, 3] => 5040,
                [1, 2, 3] => 2903040,
                [0, 0, 4] => 720,
                [0, 1, 4] => 5040,
                [1, 1, 4] => 64 * 5040,
                [0, 2, 4] => 40320,
                [1, 2, 4] => E8_SIZE,
                _ => unreachable!(),
            };
        }

        order
    }

    pub fn vertex(self) -> Point {
        let cvec: Vec8 = Mirror::ALL
            .map(|mirror| if self.has_mirror(mirror) { 1 } else { 0 })
            .into();
        let cmat = matrix![
            0, 0, 0, 0, 0, 0, -2, 2;
            0, 0, 0, 0, 0, -2, -2, 4;
            0, 0, 0, 0, -2, -2, -2, 6;
            0, 0, 0, -2, -2, -2, -2, 8;
            0, 0, 0, 0, 0, 0, 0, 4;
            -1, -1, -1, -1, -1, -1, -1, 7;
            1, -1, -1, -1, -1, -1, -1, 5;
            0, 0, -2, -2, -2, -2, -2, 10
        ];
        Point::new(cvec * cmat)
    }

    pub fn vertex_orbits(self) -> FxHashMap<Point, E8> {
        let total_vertices = E8_SIZE / self.complement().order();
        let vertex = self.vertex().representative();
        let mut seen_vertices = 0;
        let mut orbits = HashMap::from_iter([]);
        let mut rng = rand::rng();
        loop {
            let e8: E8 = rng.random();
            let v = (vertex * e8).representative();
            if orbits.insert(v, e8).is_none() {
                seen_vertices += v.orbit_size();
            }
            if seen_vertices == total_vertices {
                break;
            }
        }
        orbits
    }
}

impl Distribution<E8> for MirrorSet {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> E8 {
        let mut gens = Mirror::ALL.map(|mirror| {
            if self.has_mirror(mirror) {
                mirror.mat()
            } else {
                E8::identity()
            }
        });
        for _ in 0..100 {
            let ind1 = rng.random_range(0..8);
            let ind2 = rng.random_range(0..7);
            let ind2 = ind2 + if ind2 >= ind1 { 1 } else { 0 };
            let mul = if rng.random() {
                gens[ind2]
            } else {
                gens[ind2].inv()
            };
            if rng.random() {
                gens[ind1] = gens[ind1] * mul;
            } else {
                gens[ind1] = mul * gens[ind1];
            };
        }
        gens[rng.random_range(0..8)]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_sizes_divisible() {
        for bits_out in 0..=255 {
            for bits_in in 0..=255 {
                if bits_in & bits_out == bits_in {
                    let size_out = MirrorSet::from_bits(bits_out).unwrap().order();
                    let size_in = MirrorSet::from_bits(bits_in).unwrap().order();
                    assert_eq!(
                        size_out % size_in,
                        0,
                        "{:08b} [{}] % {:08b} [{}]",
                        bits_out,
                        size_out,
                        bits_in,
                        size_in
                    );
                }
            }
        }
    }

    #[test]
    fn all_vertex_dots_equal() {
        for bits in 1..=255 {
            let mirrors = MirrorSet::from_bits(bits).unwrap();
            let vertex = mirrors.vertex().vec();
            let dot = Mirror::ALL
                .iter()
                .find(|mirror| mirrors.has_mirror(**mirror))
                .unwrap()
                .pole()
                .dot(&vertex);
            assert_ne!(dot, 0);
            assert!(
                Mirror::ALL.iter().all(|mirror| mirror.pole().dot(&vertex)
                    == if mirrors.has_mirror(*mirror) { dot } else { 0 }),
                "{:08b}, {:?}",
                bits,
                Mirror::ALL.map(|mirror| mirror.pole().dot(&vertex)),
            )
        }
    }
}
