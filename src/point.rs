use nalgebra::RowSVector;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Point(pub RowSVector<i8, 8>);

impl Point {
    pub fn dot(self, other: Self) -> i8 {
        self.0.dot(&&other.0)
    }

    pub fn d8_orbit(self) -> D8Orbit {
        let mut coords = self.0.abs();
        coords.data.0.sort();
        D8Orbit {
            coords,
            sign: self.0.iter().fold(1, |x, y| x * y),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct D8Orbit {
    coords: RowSVector<i8, 8>,
    sign: i8,
}
