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

    pub fn components(self) -> Vec<Self> {
        let mut components = Vec::new();
        let mut middle = Self::M;

        {
            let (comps, mid) = match (self.a0(), self.a1(), self.a2(), self.a3(), self.m()) {
                (oo, oo, oo, oo, _) => (vec![], Self::empty()),
                (XX, oo, oo, oo, _) => (vec![Self::A0], Self::empty()),
                (oo, XX, oo, oo, _) => (vec![Self::A1], Self::empty()),
                (XX, XX, oo, oo, _) => (vec![Self::A0 | Self::A1], Self::empty()),
                (oo, oo, XX, oo, _) => (vec![Self::A2], Self::empty()),
                (XX, oo, XX, oo, _) => (vec![Self::A0, Self::A2], Self::empty()),
                (oo, XX, XX, oo, _) => (vec![Self::A1 | Self::A2], Self::empty()),
                (XX, XX, XX, oo, _) => (vec![Self::A0 | Self::A1 | Self::A2], Self::empty()),
                (oo, oo, oo, XX, oo) => (vec![Self::A3], Self::empty()),
                (XX, oo, oo, XX, oo) => (vec![Self::A0, Self::A3], Self::empty()),
                (oo, XX, oo, XX, oo) => (vec![Self::A1, Self::A3], Self::empty()),
                (XX, XX, oo, XX, oo) => (vec![Self::A0 | Self::A1, Self::A3], Self::empty()),
                (oo, oo, XX, XX, oo) => (vec![Self::A2 | Self::A3], Self::empty()),
                (XX, oo, XX, XX, oo) => (vec![Self::A0, Self::A2 | Self::A3], Self::empty()),
                (oo, XX, XX, XX, oo) => (vec![Self::A1 | Self::A2 | Self::A3], Self::empty()),
                (XX, XX, XX, XX, oo) => (
                    vec![Self::A0 | Self::A1 | Self::A2 | Self::A3],
                    Self::empty(),
                ),
                (oo, oo, oo, XX, XX) => (vec![], Self::A3),
                (XX, oo, oo, XX, XX) => (vec![Self::A0], Self::A3),
                (oo, XX, oo, XX, XX) => (vec![Self::A1], Self::A3),
                (XX, XX, oo, XX, XX) => (vec![Self::A0 | Self::A1], Self::A3),
                (oo, oo, XX, XX, XX) => (vec![], Self::A2 | Self::A3),
                (XX, oo, XX, XX, XX) => (vec![Self::A0], Self::A2 | Self::A3),
                (oo, XX, XX, XX, XX) => (vec![], Self::A1 | Self::A2 | Self::A3),
                (XX, XX, XX, XX, XX) => (vec![], Self::A0 | Self::A1 | Self::A2 | Self::A3),
            };
            components.extend(comps);
            middle |= mid;
        }

        {
            let (comps, mid) = match (self.b0(), self.b1(), self.m()) {
                (oo, oo, oo) => (vec![], Self::empty()),
                (XX, oo, oo) => (vec![Self::B0], Self::empty()),
                (oo, XX, oo) => (vec![Self::B1], Self::empty()),
                (XX, XX, oo) => (vec![Self::B0 | Self::B1], Self::empty()),
                (oo, oo, XX) => (vec![], Self::empty()),
                (XX, oo, XX) => (vec![Self::B0], Self::empty()),
                (oo, XX, XX) => (vec![], Self::B1),
                (XX, XX, XX) => (vec![], Self::B0 | Self::B1),
            };
            components.extend(comps);
            middle |= mid;
        }

        {
            let (comps, mid) = match (self.c(), self.m()) {
                (oo, oo) => (vec![], Self::empty()),
                (XX, oo) => (vec![Self::C], Self::empty()),
                (oo, XX) => (vec![], Self::empty()),
                (XX, XX) => (vec![], Self::C),
            };
            components.extend(comps);
            middle |= mid;
        }

        if self.m() == XX {
            components.push(middle);
        }

        components
    }

    fn component_order(self) -> u64 {
        if self.a3() == XX && self.b1() == XX && self.c() == XX {
            // branched component
            if self.b0() == XX {
                // E component
                if self.a0() == XX {
                    696729600
                } else if self.a1() == XX {
                    2903040
                } else if self.a2() == XX {
                    51840
                } else {
                    1920
                }
            } else {
                // D component
                if self.a0() == XX {
                    322560
                } else if self.a1() == XX {
                    23040
                } else if self.a2() == XX {
                    1920
                } else {
                    192
                }
            }
        } else {
            // A component
            match self.0.0.count_ones() {
                1 => 2,
                2 => 6,
                3 => 24,
                4 => 120,
                5 => 720,
                6 => 5040,
                _ => 40320, // 7
            }
        }
    }

    pub fn order(self) -> u64 {
        self.components()
            .into_iter()
            .map(Self::component_order)
            .product()
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
