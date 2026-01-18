use nalgebra::RowSVector;
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
        Self(self.0.map(|p| {
            let other_p = other.0[p.ax()];
            AxSign::new(other_p.ax(), p.sign() * other_p.sign())
        }))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Point {
    pub sorted: Vec8,
    pub sign: i16,
    pub d8: D8,
}

impl Point {
    pub fn new(vec: Vec8) -> Self {
        let sign = vec.iter().map(|x| x.signum()).product();
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
        deco.sort();
        let sorted = deco.map(|(x, _)| x).into();
        let mut d8 = D8(deco.map(|(_, p)| p));
        d8.0[0].flip_sign();
        Self {
            sorted,
            sign,
            d8: d8.inv(),
        }
    }

    pub fn vec(self) -> Vec8 {
        let mut rep = self.sorted;
        rep[0] *= self.sign;
        rep * self.d8
    }

    pub fn dot(self, other: Self) -> i16 {
        self.vec().dot(&other.vec())
    }

    // pub fn d8_orbit(self) -> D8Orbit {
    //     let mut coords = self.0.abs();
    //     coords.data.0.sort();
    //     D8Orbit {
    //         coords,
    //         sign: self.0.iter().fold(1, |x, y| x * y),
    //     }
    // }
}

impl Mul<D8> for Point {
    type Output = Self;
    fn mul(self, d8: D8) -> Self {
        Self {
            sorted: self.sorted,
            sign: self.sign,
            d8: self.d8 * d8,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn d8_order() {
        let d8 = D8([1, 2, !0, 3, 4, 5, 6, !7].map(AxSign));
        assert_eq!(d8 * d8 * d8 * d8 * d8 * d8, D8::identity());
        let d8_2 = D8([0, 1, 2, !6, 3, 4, 5, !7].map(AxSign));
        assert_eq!(d8 * d8.inv(), D8::identity());
        assert_eq!(d8_2 * d8_2.inv(), D8::identity());
        assert_eq!(d8 * d8_2 * d8.inv() * d8_2.inv(), D8::identity());
    }

    #[test]
    fn d8_roundtrip() {
        let v: Vec8 = [10, 4, 6, 2, -2, 8, 6, 2].into();
        let point = Point::new(v);
        assert_eq!(point.vec(), v, "{:?}", point);
    }
}
