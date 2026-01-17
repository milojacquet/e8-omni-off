use nalgebra::RowSVector;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Point(pub RowSVector<i8, 8>);

impl Point {
    pub fn dot(self, other: Self) -> i8 {
        self.0.dot(&&other.0)
    }
}
