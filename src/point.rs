use crate::combs::COMBS;
use nalgebra::RowSVector;
use std::iter::once;
use std::ops::Mul;

pub type Vec8 = RowSVector<i16, 8>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AxSign(i8);

impl AxSign {
    pub fn new(ax: usize, sign: i16) -> Self {
        if sign >= 0 {
            Self(ax as i8)
        } else {
            Self(!(ax as i8))
        }
    }

    pub fn ax(self) -> usize {
        if self.0 >= 0 {
            self.0 as usize
        } else {
            (!self.0) as usize
        }
    }

    pub fn sign(self) -> i16 {
        if self.0 >= 0 { 1 } else { -1 }
    }

    pub fn flip_sign(&mut self) {
        self.0 = !self.0;
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct D8([AxSign; 8]);

impl D8 {
    pub fn new(axs: [AxSign; 8]) -> Self {
        let new = Self(axs);
        assert!(new.signs_even(), "{:?}", axs);
        new
    }

    pub fn identity() -> Self {
        Self([0, 1, 2, 3, 4, 5, 6, 7].map(AxSign))
    }

    pub fn inv(self) -> Self {
        let mut inv = [AxSign(-128); 8];
        for (i, p) in self.0.iter().enumerate() {
            inv[p.ax()] = AxSign::new(i, p.sign());
        }
        Self(inv)
    }

    pub fn signs_even(self) -> bool {
        self.0.iter().filter(|x| x.sign() == -1).count() % 2 == 0
    }
}

impl Mul<D8> for Vec8 {
    type Output = Self;
    fn mul(self, d8: D8) -> Self {
        d8.0.map(|p| self[p.ax()] * p.sign()).into()
    }
}

impl Mul for D8 {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self(other.0.map(|p| {
            let self_p = self.0[p.ax()];
            AxSign::new(self_p.ax(), p.sign() * self_p.sign())
        }))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Orbit {
    pub rep: Vec8,
    pub sign: i16,
}

impl Orbit {
    pub fn size(self) -> u64 {
        let mut size = 1;
        let mut group = 1;
        for (i, (&x, y)) in self
            .rep
            .iter()
            .zip(self.rep.iter().skip(1).map(Some).chain(once(None)))
            .enumerate()
        {
            if Some(&x) == y {
                group += 1;
            } else {
                size *= COMBS[i + 1][group].len() as u64;
                group = 1;
            }
            if x != 0 && i != 0 {
                size *= 2;
            }
        }
        size
    }

    pub fn iter(self) -> impl Iterator<Item = Point> {
        let mut bases = Vec::new();
        let mut group = 1;
        let mut signs = 0;
        for (i, (&x, y)) in self
            .rep
            .iter()
            .zip(self.rep.iter().skip(1).map(Some).chain(once(None)))
            .enumerate()
        {
            if Some(&x) == y {
                group += 1;
            } else {
                bases.push((i + 1, group, COMBS[i + 1][group].len() as u64));
                group = 1;
            }
            if x != 0 && i != 0 {
                signs += 1;
            }
        }
        bases.reverse();
        let bases = bases;

        let perms = bases.iter().map(|(_, _, p)| p).product();

        (0..perms)
            .map(move |mut i| {
                let mut d8 = [None; 8];
                for (n, k, p) in &bases {
                    let comb_num = COMBS[*n][*k][(i % p) as usize];
                    i /= p;

                    let none_iter = d8.iter_mut().filter(|p| p.is_none());
                    for (j, (_, p)) in none_iter
                        .enumerate()
                        .filter(|(b, _)| (comb_num >> (n - 1 - b)) & 1 == 1)
                        .enumerate()
                    {
                        *p = Some(AxSign::new(n - k + j, 1)) // placeholder sign
                    }
                }

                let d8 = d8.map(Option::unwrap);

                (0u8..1 << signs).map(move |s| {
                    let mut d8 = d8;
                    for p in d8.iter_mut() {
                        if (s >> (7 - p.ax())) & 1 == 1 || (p.ax() == 0 && s.count_ones() % 2 != 0)
                        {
                            p.flip_sign();
                        }
                    }
                    Point {
                        orbit: self,
                        d8: D8::new(d8),
                    }
                })
            })
            .flatten()
    }
}

pub fn opt_bits_to_num(bits: [Option<bool>; 8]) -> u8 {
    let mut num = 0;
    for bit in bits {
        if let Some(bit) = bit {
            num <<= 1;
            num |= if bit { 1 } else { 0 };
        }
    }
    num
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Point {
    pub orbit: Orbit,
    pub d8: D8,
}

impl Point {
    pub fn new(vec: Vec8) -> Self {
        let sign = vec.iter().map(|x| x.signum()).product();
        let sign_pm1: i16 = vec
            .iter()
            .map(|x| if x.signum() == 0 { 1 } else { x.signum() })
            .product();
        let mut deco = [
            (vec[0].abs(), AxSign::new(0, vec[0].signum())),
            (vec[1].abs(), AxSign::new(1, vec[1].signum())),
            (vec[2].abs(), AxSign::new(2, vec[2].signum())),
            (vec[3].abs(), AxSign::new(3, vec[3].signum())),
            (vec[4].abs(), AxSign::new(4, vec[4].signum())),
            (vec[5].abs(), AxSign::new(5, vec[5].signum())),
            (vec[6].abs(), AxSign::new(6, vec[6].signum())),
            (vec[7].abs(), AxSign::new(7, vec[7].signum())),
        ]; // why won't they stabilize zip
        deco.sort_by_key(|(x, _)| *x);
        let rep = deco.map(|(x, _)| x).into();
        let mut d8_inner = deco.map(|(_, p)| p);
        if sign_pm1 == -1 {
            d8_inner[0].flip_sign();
        }
        // eprintln!("{:?}", vec);
        let d8 = D8::new(d8_inner);
        Self {
            orbit: Orbit { rep, sign },
            d8: d8.inv(),
        }
    }

    pub fn vec(self) -> Vec8 {
        let mut rep = self.orbit.rep;
        rep[0] *= self.orbit.sign;
        rep * self.d8
    }

    pub fn dot(self, other: Self) -> i16 {
        self.vec().dot(&other.vec())
    }

    pub fn orbit_size(self) -> u64 {
        self.orbit.size()
    }

    pub fn orbit_index(self) -> u64 {
        let mut index = 0;

        let mut group = 1;
        let mut signs = 0;
        for (i, (&x, y)) in self
            .orbit
            .rep
            .iter()
            .zip(self.orbit.rep.iter().skip(1).map(Some).chain(once(None)))
            .enumerate()
        {
            if Some(&x) == y {
                group += 1;
            } else {
                let comb_num = opt_bits_to_num(self.d8.0.map(|p| {
                    if p.ax() > i {
                        None
                    } else {
                        Some(p.ax() + group > i)
                    }
                }));
                index *= COMBS[i + 1][group].len() as u64;
                index += COMBS[i + 1][group].binary_search(&comb_num).unwrap() as u64;
                group = 1;
            }
            if x != 0 && i != 0 {
                signs += 1;
            }
        }

        index <<= signs;
        for p in self.d8.0 {
            if p.ax() != 0 && self.orbit.rep[p.ax()] != 0 && p.sign() == -1 {
                index |= 1 << (7 - p.ax())
            }
        }

        index
    }

    pub fn representative(self) -> Self {
        Self {
            orbit: self.orbit,
            d8: D8::identity(),
        }
    }
}

impl Mul<D8> for Point {
    type Output = Self;
    fn mul(self, d8: D8) -> Self {
        // TODO: don't roundtrip through vec
        Self::new(
            Self {
                orbit: self.orbit,
                d8: self.d8 * d8,
            }
            .vec(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn d8_order() {
        let d8 = D8::new([1, 2, !0, 3, 4, 5, 6, !7].map(AxSign));
        assert_eq!(d8 * d8 * d8 * d8 * d8 * d8, D8::identity());
        let d8_2 = D8::new([0, 1, 2, !6, 3, 4, 5, !7].map(AxSign));
        assert_eq!(d8 * d8.inv(), D8::identity());
        assert_eq!(d8_2 * d8_2.inv(), D8::identity());
        assert_eq!(d8 * d8_2 * d8.inv() * d8_2.inv(), D8::identity());
    }

    #[test]
    fn d8_roundtrip_1() {
        let v: Vec8 = [10, 4, 6, 2, -2, 8, 6, 2].into();
        let point = Point::new(v);
        assert_eq!(point.vec(), v, "{:?}", point);
    }

    #[test]
    fn d8_roundtrip_2() {
        let v: Vec8 = [1, -1, -1, -1, -1, -1, -1, 5].into();
        let point = Point::new(v);
        assert_eq!(point.vec(), v, "{:?}", point);
    }

    #[test]
    fn d8_orbit_size() {
        assert_eq!(
            Orbit {
                rep: [1, 2, 3, 4, 5, 6, 7, 8].into(),
                sign: 1
            }
            .size(),
            128 * 40320
        );
        assert_eq!(
            Orbit {
                rep: [0, 2, 3, 4, 5, 6, 7, 8].into(),
                sign: 1
            }
            .size(),
            128 * 40320
        );
        assert_eq!(
            Orbit {
                rep: [0, 0, 3, 4, 5, 6, 7, 8].into(),
                sign: 1
            }
            .size(),
            64 * 40320 / 2
        );
        assert_eq!(
            Orbit {
                rep: [1, 1, 1, 1, 1, 1, 1, 1].into(),
                sign: 1
            }
            .size(),
            128
        );
    }

    fn d8_orbit_index_consistent(rep: [i16; 8]) {
        let orbit = Orbit {
            rep: rep.into(),
            sign: 1,
        };
        for (i, point) in orbit.iter().enumerate() {
            if i % 100000 == 0 {
                println!("{i}")
            }
            assert_eq!(point.orbit_index(), i as u64, "{}", point.vec());
        }
    }

    #[ignore = "very long"]
    #[test]
    fn d8_orbit_index_consistent_12345678() {
        d8_orbit_index_consistent([1, 2, 3, 4, 5, 6, 7, 8])
    }

    #[ignore = "very long"]
    #[test]
    fn d8_orbit_index_consistent_01234567() {
        d8_orbit_index_consistent([0, 1, 2, 3, 4, 5, 6, 7])
    }

    #[test]
    fn d8_orbit_index_consistent_12233344() {
        d8_orbit_index_consistent([1, 2, 2, 3, 3, 3, 4, 4])
    }

    #[test]
    fn d8_orbit_index_consistent_00122333() {
        d8_orbit_index_consistent([0, 0, 1, 2, 2, 3, 3, 3])
    }

    #[test]
    fn d8_mul_weird() {
        let v: Vec8 = [0, 0, 0, 0, 0, 0, 2, 2].into();
        let d8_1 = D8::new([
            AxSign(0),
            AxSign(1),
            AxSign(2),
            AxSign(3),
            AxSign(-7),
            AxSign(4),
            AxSign(5),
            AxSign(-8),
        ]);
        let d8_2 = D8::new([
            AxSign(-2),
            AxSign(2),
            AxSign(-5),
            AxSign(-6),
            AxSign(-7),
            AxSign(0),
            AxSign(-4),
            AxSign(-8),
        ]);

        assert_eq!(v * (d8_1 * d8_2), (v * d8_1) * d8_2)
    }

    #[test]
    fn point_mul_eq() {
        let point = Point::new([0, 0, 0, 0, 0, 0, 2, 2].into());
        let d8 = D8::new([1, 0, 2, 3, 4, 5, 6, 7].map(AxSign));
        assert_eq!(point * d8, point);
    }

    #[test]
    fn point_d8_even() {
        let point = Point::new([0, 0, 0, 0, 0, 0, -2, 2].into());
        assert!(point.d8.signs_even());
    }
}
